# StochasticOptimalControl_ToyExample

NOTE: Running the code requires installation of CasADi (https://web.casadi.org/).

Here we present a stochastic optimal control approach that allows performing simulations in the presence of noise and state uncertainty with non-linear systems.
To bring the focus to the implementation as much as possible, the selected example has linear dynamics, however non-linear constraints path constraints are presented.

The simulated movement can be thought of as a point to point reaching movement whereby an obstacle needs to be avoided. A point-mass is actuated in a 2D plane by a horizontal and vertical force:

- In the deterministic case the energetically optimal controller will take the shortest path and 'chafe' the obstacle. The controller here is open-loop.

- In the presence of motor noise and an open-loop controller avoiding the obstacle with a specified certainty (e.g. 95%) requires offsetting the trajectory from the obstacle.

- In the presence of motor noise and a feedback controller avoiding the obstacle with a specified certainty (e.g. 95%) requires offsetting the trajectory from the obstacle. The offset can be less however than in the open-loop controller. Under a feedback controller, the end-point variability can be controlled as well.


Optimal open-loop deterministic trajectory and forward simulations with optimal deterministic controller in a stochastic environment:
![alt text](https://github.com/tomvanwouwe1992/StochasticOptimalControl_ToyExample/blob/master/nominalSolution.jpg)


Optimal open-loop and feedback controlled stochastic trajectory and forward simulations with optimal controllers in a stochastic environment:
![alt text](https://github.com/tomvanwouwe1992/StochasticOptimalControl_ToyExample/blob/master/stochasticSolutionOLvsFB.jpg)
