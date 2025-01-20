The magnetic mirror corresponds to  ... <br>

Magnetic trap using a Bottle profile




https://github.com/user-attachments/assets/38e82c09-b26f-4e67-afc2-482a41cbcfb0





 Let's start with a simpler model, where a particle with mass m and charge q is under the effect of a magentic field where the z profile is given by

$$B_z=B_{z0}(1+z)$$

...

When the magnetic field becomes strong enough, the parallel velocity ($v_\parallel$) of the particle goes to zero, and the perpendicular velocity ($v_\perp$) reaches its maximum value. At this point, the particle is reflected and begins moving in the opposite direction along the  z-axis (parallel to the magnetic field lines). 






https://github.com/user-attachments/assets/a51226dc-3823-4d79-b7e9-62bacf5c6a37






In a very slowly magnetic field, the magnetic moment $\mu=mv^2_\perp/2B$ is an adiabatic constant, so it leaves the relation <br>

$$v^2_\perp \propto 2B$$

such that for regions where the magnetic field increases, it also does the perpendicular velocity. Additionally, the energy E is also conserved and is purelly kinetic given by 

$$ E= \frac{1}{2} m ( v^2_\parallel  + v^2_\perp)$$



https://github.com/user-attachments/assets/4d83f835-c17f-469f-8cc9-edb52791b311

When $v^2_\perp$ increases, $v^2_\parallel$ must decrease to conserve the particle's total kinetic energy. Reflection occurs when 
$v^2_\parallel=0$, meaning the particle momentarily comes to a stop along the z-axis. At this point, all its kinetic energy is concentrated in the perpendicular direction ($v^2_\perp$ is maximized). Afterward, the particle reverses its direction of motion along the z-axis

![height_velocities](https://github.com/user-attachments/assets/2b96554e-c717-4b1a-beee-63ce745fa2e9)

 The gradient of the magnetic field ($\nabla B$) plays a crucial role in determining the system's behavior. It governs the adiabaticity of a particle's motion and the conservation of its magnetic moment. Naturally, this raises the question: how small must the gradient be to ensure reflection? The answer is that the magnetic field must vary slowly relative to the particle's motion and we can quantify this behavior by the adiabatic parameter $\epsilon$ defined as

$$\epsilon = \frac{r_L}{L} <1$$

where $r_L=m v_\parallel/qB$ is the Larmor radius and $L=B/\nabla B$ is the magnetic field gradient scale length.
 
<img width="448" alt="epsilon" src="https://github.com/user-attachments/assets/96536ab7-3ae0-4363-b796-1f474c216b01" />


Symplectic integrator:
1. Leap frog...
2. Verlet...
