# Motion Planning in Frenet-Serret Frame
This repo provides modular codebase in both python and cpp, to track a mule-led-trailer robot in any path of random shape, given its parametric equation. The path is converted to the Frenet-Serret frame and lateral deviation from reference path is penalized, to achieve strict adherance. A Model Predictive Control strategy is implemented to for trajectory tracking. The robot dynamics and simulation outputs are provided in further sections. 

## Requirements
* Python 3.13.5
* matplotlib 3.10.3
* scipy 1.6.1

## Robot Dynamics
The position of the robot (mule) is given by $(x, y, θ)$. This point is at the centre of the front axle. Along the line from this point to the rear wheel and beyond is the hitch, at distance $L$. Further, we have a trailer connected to this hitch point. The rear wheel of this trailer is at a distance $D$ from the hitch. We assume the line connecting the hitch and the two back wheels form an angle $\phi$.
The kinematics equations of the robot is the well-known unicycle model:<br>
>$\dot x = v cos \theta$<br>
>$\dot y = v sin \theta$<br>
>$\dot \theta = \omega$<br>

In turn, the kinematics equations of the trailer can be derived from its linear ($v_{t}$) and angular velocities ($\dot\phi$). This is given by:<br>
>$v_{t} = v cos (\theta-\phi) + L\omega sin (\theta-\phi)$<br>
>$\dot\phi = v sin (\theta-\phi) - L\omega cos (\theta-\phi)$

## Simulations
* Path Following on Figure Eight Curvatue
![](path_following.gif)

* Obstacle Avoidance on Figure Eight Curvatue
![](obstacle_avoidance.gif)

* Reverse Path Tracking on Figure Eight Curvatue
![](reverse_tracking.gif)

## Play Around
* Try other inbuilt parametric paths like circle, ellipse or even a straight line!
* Increasing $N$ gives better results in terms of accuracy.
* Weights can be tuned further for stricter path following.

## NOTE
* Obstacle avoidance is not yet supported for reverse tracking motion.

