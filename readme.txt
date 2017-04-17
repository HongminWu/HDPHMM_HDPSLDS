@ Hongmin Wu
1. This is based on software originally written by Emily B. Fox and adapted for using in the fault detection and recognition in robot contract task.
    if you want to get further informations, please don't hesitate to contact me in Github/Email.
2. The main script to run is "demos" folder, there are two phase to accomplish the fault detection and recognition, 
   first, demos/demo_training: the script denotes that training the specific state by nonparametric hmm and inference the unknown parameters
   second, demos/demo_testing: we use the log-likelihood of observations to calssify the success/failure by calculate the p(y_{i}|trained_HMM_model)
   finally, global_variables : the global parameters should be defined before run any demos, if you want to create a struct that has new observations
   of each of your tasks, you can following this file to add appropriate parameters.

3. Why we  are crazy about this research?
this work investigates the use of Bayesian nonparametric techniques that define a prior on Markov Switching processes with an unbounded number of states for the state estimation of contact tasks. In particular, we explore the use of a sticky Hierarchical Dirichlet stochastic process prior to learn a number of Markov Jump Linear System models for wrench signatures in the contact tasks. We believe the approach will lead to strong
generalizations that discover and model the underlying robot state in robust, accurate,
and computationally efficient ways. Previous state estimation techniques for contact
tasks have often used Markov Models or supervised learning methods. Both of these methods
have primarily been parametric, implying that the number of classes or states is defined
a priori and fixed through the entire task. This is a limiting assumption as it lacks
modelling power to properly segment the data into natural clusters that represent task’s
state. Furthermore, Markov Models which are a type of Markov Switching Process, are
limited in their ability to model complex dynamical behaviors due to their assumption
that observations are conditionally independent given the state. Instead, this work
considers more powerful modelling Markov Switching Processes; namely, those in the class
of Markov Jump Linear Systems (MJLS). These systems more powerfully describe rich
dynamical phenomena through their conditionally linear modes. In particular, the
switching vector autoregressive and the switching linear dynamical system methods will be
used. Key scientific issues to be solved lie in the determination of dynamic parameters
for the MJLS methods and the determination of more informative prior distributions and
their corresponding hyper-parameters. These methods will be tested on contact tasks that
consist of multiple contact sub-tasks including pushing, pulling, alignments, insertions,
mating, sliding, amongst others. We expect to integrate a real-time open source framework
that will exploit MJLS system and significantly increase the robustness, accuracy, and
computational efficiency of state-of-the art approaches in robot contact task state
estimation. We also expect to compile and open to the public a library of robot primitive
actions useful in contact state estimation. A wide array of robotic domains will benefit
from the ability to confirm commanded actions including: collaborative manufacturing
robots, service robots, as well as areas like high-level planning and programming by
demonstration. Better state estimates will also derive in error reduction and increased
efficiency across all robot domains.
4. we have some datesets of HIRO robot snap assembly can be download at: http://www.juanrojas.net/
