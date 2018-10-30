
import argparse

import Tkinter as tk
import numpy as np

from pydrake.examples.manipulation_station import ManipulationStation
from pydrake.geometry import ConnectDrakeVisualizer
from pydrake.multibody.multibody_tree.multibody_plant import MultibodyPlant
from pydrake.manipulation.simple_ui import SchunkWsgButtons
from pydrake.manipulation.planner import (
    DifferentialInverseKinematicsParameters, DoDifferentialInverseKinematics)
from pydrake.math import RigidTransform, RollPitchYaw
from pydrake.systems.analysis import Simulator
from pydrake.systems.framework import ( AbstractValue, BasicVector,
                                        DiagramBuilder,
                                        LeafSystem,
                                        PortDataType )
from pydrake.util.eigen_geometry import Isometry3, AngleAxis


class EndEffectorTeleop(LeafSystem):
    def __init__(self):
        LeafSystem.__init__(self)
        self._DeclareAbstractOutputPort("X_WE",
              lambda: AbstractValue.Make(Isometry3.Identity()),
              self._DoCalcOutput)

        self._DeclarePeriodicPublish(0.1, 0.0)

        self.window = tk.Tk()
        self.window.title("End-Effector TeleOp")

        self.roll  = tk.Scale(self.window, from_=-2 * np.pi, to=2 * np.pi,
                              resolution=-1,
                              label="roll",
                              length=800,
                              orient=tk.HORIZONTAL)
        self.roll.pack()
        self.pitch = tk.Scale(self.window, from_=-2 * np.pi, to=2 * np.pi,
                              resolution=-1,
                              label="pitch",
                              length=800,
                              orient=tk.HORIZONTAL)
        self.pitch.pack()
        self.yaw   = tk.Scale(self.window, from_=-2 * np.pi, to=2 * np.pi,
                              resolution=-1,
                              label="yaw",
                              length=800,
                              orient=tk.HORIZONTAL)
        self.yaw.pack()
        self.x     = tk.Scale(self.window, from_=0.1, to=0.8,
                              resolution=-1,
                              label="x",
                              length=800,
                              orient=tk.HORIZONTAL)
        self.x.pack()
        self.y     = tk.Scale(self.window, from_=-0.3, to=0.3,
                              resolution=-1,
                              label="y",
                              length=800,
                              orient=tk.HORIZONTAL)
        self.y.pack()
        self.z     = tk.Scale(self.window, from_=0, to=0.8,
                              resolution=-1,
                              label="z",
                              length=800,
                              orient=tk.HORIZONTAL)
        self.z.pack()

    # @param pose is an Isometry3
    def SetPose(self, pose):
        tf = RigidTransform(pose)
        self.SetRPY(RollPitchYaw(tf.rotation()))
        self.SetXYZ(pose.translation())

    # @param rpy is a 3 element vector of roll, pitch, yaw
    def SetRPY(self, rpy):
        self.roll.set(rpy.roll_angle())
        self.pitch.set(rpy.pitch_angle())
        self.yaw.set(rpy.yaw_angle())

    def SetXYZ(self, xyz):
        self.x.set(xyz[0])
        self.y.set(xyz[1])
        self.z.set(xyz[2])

    def _DoPublish(self, context, event):
        self.window.update_idletasks()
        self.window.update()

    def _DoCalcOutput(self, context, output):
        output.get_mutable_value().set_rotation(RollPitchYaw(self.roll.get(),
                                                             self.pitch.get(),
                                                             self.yaw.get()).
                                                ToRotationMatrix().matrix())
        output.get_mutable_value().set_translation([self.x.get(), self.y.get(),
                                                    self.z.get()])

class DifferentialIK(LeafSystem):
    # @param robot is a reference to a MultibodyPlant.
    # @param frame_E is a multibody::Frame on the robot.
    # @param params is a DifferentialIKParams.
    def __init__(self, robot, frame_E, parameters, time_step):
        LeafSystem.__init__(self)
        self.robot = robot
        self.frame_E = frame_E
        self.parameters = parameters
        self.parameters.set_timestep(time_step)
        self.time_step = time_step
        self.robot_context = robot.CreateDefaultContext()
        x = self.robot.tree().get_mutable_multibody_state_vector(
            self.robot_context)
        x[:] = 0

        # Store the robot positions as state.
        self._DeclareDiscreteState(robot.num_positions())
        self._DeclarePeriodicDiscreteUpdate(time_step)

        # Desired pose of frame E in world frame.
        self._DeclareAbstractInputPort("X_WE_desired")

        # Alternatively, provide the output as desired positions.
        self._DeclareVectorOutputPort("desired_position", BasicVector(
            robot.num_positions()), self.CopyPositionOut)

    def SetPositions(self, context, q):
        context.get_mutable_discrete_state(0).SetFromVector(q)

    def ForwardKinematics(self, q):
        x = self.robot.tree().get_mutable_multibody_state_vector(
            self.robot_context)
        x[:robot.num_positions()] = q
        return self.robot.tree().EvalBodyPoseInWorld(
            self.robot_context, self.frame_E.body())

    def CalcPoseError(self, X_WE_desired, q):
        pose = self.ForwardKinematics(q)
        err_vec = np.zeros(6)
        err_vec[-3:] = X_WE_desired.translation() - pose.translation()

        rot_err = AngleAxis(X_WE_desired.rotation() * pose.rotation().transpose())
        err_vec[:3] = rot_err.axis() * rot_err.angle()

    def _DoCalcDiscreteVariableUpdates(
            self, context, events, discrete_state):
        X_WE_desired = self.EvalAbstractInput(context, 0).get_value()
        q_last = context.get_discrete_state_vector().get_value()

        x = self.robot.tree().get_mutable_multibody_state_vector(
            self.robot_context)
        x[:robot.num_positions()] = q_last
        result = DoDifferentialInverseKinematics(self.robot,
                                                self.robot_context,
                                                X_WE_desired, self.frame_E,
                                                self.parameters)

        if (result.status != result.status.kSolutionFound):
            discrete_state.get_mutable_vector().SetFromVector(q_last)
        else:
            discrete_state.get_mutable_vector().SetFromVector(q_last +
                                   self.time_step*result.joint_velocities)

    def CopyPositionOut(self, context, output):
        output.SetFromVector(context.get_discrete_state_vector().get_value())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--target_realtime_rate", type=float, default=1.0,
    help="Desired rate relative to real time.  See documentation for "
         "Simulator::set_target_realtime_rate() for details.")
parser.add_argument(
    "--duration", type=float, default=np.inf,
    help="Desired duration of the simulation in seconds.")
args = parser.parse_args()

builder = DiagramBuilder()

time_step = 0.002
station = builder.AddSystem(ManipulationStation(time_step))
station.AddCupboard()
station.Finalize()
q0 = [0, 0.6, 0, -1.75, 0, 1.0, 0]

robot = station.get_controller_plant()
params = DifferentialInverseKinematicsParameters(robot.num_positions(),
                                                 robot.num_velocities())

params.set_timestep(time_step)
params.set_nominal_joint_position(q0)
iiwa14_velocity_limits = np.array([1.4, 1.4, 1.7, 1.3, 2.2, 2.3, 2.3])
params.set_joint_velocity_limits((-.5*iiwa14_velocity_limits,
                                  .5*iiwa14_velocity_limits))

differential_ik = builder.AddSystem(DifferentialIK(
    robot, robot.GetFrameByName("iiwa_link_7"), params, time_step))

builder.Connect(differential_ik.GetOutputPort("desired_position"),
 station.GetInputPort("iiwa_position"))

teleop = builder.AddSystem(EndEffectorTeleop())
builder.Connect(teleop.get_output_port(0),
                differential_ik.GetInputPort("X_WE_desired"))

#wsg_buttons = builder.AddSystem(SchunkWsgButtons(teleop.window))
#builder.Connect(wsg_buttons.GetOutputPort("position"), station.GetInputPort(
#    "wsg_position"))
#builder.Connect(wsg_buttons.GetOutputPort("force_limit"),
#                station.GetInputPort("wsg_force_limit"))

ConnectDrakeVisualizer(builder, station.get_mutable_scene_graph(),
                       station.GetOutputPort("pose_bundle"))

diagram = builder.Build()
simulator = Simulator(diagram)

station_context = diagram.GetMutableSubsystemContext(station,
                                             simulator.get_mutable_context())

station.SetIiwaPosition(q0, station_context)
station.SetIiwaVelocity(np.zeros(7), station_context)
#station.SetWsgPosition(0.1, station_context)
#station.SetWsgVelocity(0, station_context)

teleop.SetPose(differential_ik.ForwardKinematics(q0))
differential_ik.SetPositions(diagram.GetMutableSubsystemContext(
    differential_ik, simulator.get_mutable_context()), q0)

station_context.FixInputPort(station.GetInputPort(
    "iiwa_feedforward_torque").get_index(), np.zeros(7))

simulator.set_target_realtime_rate(args.target_realtime_rate)
simulator.StepTo(args.duration)
