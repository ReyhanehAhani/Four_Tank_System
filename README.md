
# Controlling Four Tank System

This project is consist of three phases. In each phase we test and simulate the system using MATLAB Simulink.



 -  ## Phase 1 
  
The main objective of this phase is identifying the system, so at the first step we linearized equations of system and find the valid range for linear equations, then we check controllability and observability of system and finally we used kalman decomposition to obtain subsystems.
here we can see  linearized state space matrixes.
![Screen Shot 1402-02-21 at 22 47 22](https://github.com/ReyhanehAhani/Four_Tank_System/assets/88882191/7a43e0b3-f37e-4885-88a1-21aeaffac5f6)




- ## Phase 2
After indentifying the system and obtaining the minimal subsytem, we use state feedback controller and simulated it in simulink and then we added static pre-compensator and dynamic pre-compensator in order to make the system performes better in reference signal tracking.
Here are some pictures of adding controller and pre-compensator.

<img width="473" alt="Screen Shot 1402-02-21 at 23 17 41" src="https://github.com/ReyhanehAhani/Four_Tank_System/assets/88882191/6d94af86-609a-40c9-8334-33ece81ad90e">
<img width="610" alt="Screen Shot 1402-02-21 at 23 18 56" src="https://github.com/ReyhanehAhani/Four_Tank_System/assets/88882191/205f99c5-369b-44b2-8d1e-749d7606114e">
<img width="611" alt="Screen Shot 1402-02-21 at 23 19 38" src="https://github.com/ReyhanehAhani/Four_Tank_System/assets/88882191/e555abfa-e80a-415e-9645-ee7d7e7df4d1">
<img width="609" alt="Screen Shot 1402-02-21 at 23 20 15" src="https://github.com/ReyhanehAhani/Four_Tank_System/assets/88882191/209ab825-b525-475f-b32c-29a3b86d5cc8">
<img width="609" alt="Screen Shot 1402-02-21 at 23 20 15" src="https://github.com/ReyhanehAhani/Four_Tank_System/assets/88882191/3c1a30d9-6c26-4aad-b7f4-3e8814aaec4f">

 
 
 - ## Phase 3
 In this phase we focouse on observability and stability of this system and tried to design a full order and reduced order observer , then we check all type of stabilities and finally we designed the linear optimal control system and simulated this on non-linear system.
 Here is some pictures of observer block diagram on simulink.

<img width="913" alt="Screen Shot 1402-02-21 at 23 21 43" src="https://github.com/ReyhanehAhani/Four_Tank_System/assets/88882191/6e4ed40c-01ab-43da-9d8d-bb261eb9a882">
<img width="923" alt="Screen Shot 1402-02-21 at 23 22 38" src="https://github.com/ReyhanehAhani/Four_Tank_System/assets/88882191/f8686b0b-5654-4265-9bce-9d0d2e7e3cbb">








