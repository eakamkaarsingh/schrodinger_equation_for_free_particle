// schrodinger equation for free particle
// $\frac{-\hbar^2}{2m} \frac{\partial ^2\psi}{\partial x^2} + V\psi = i \hbar \frac{\partial \psi}{\partial t}$


// since we are talking about free particle, the potential will be zero; so we get:
// $\frac{-\hbar^2}{2m} \frac{\partial ^2\psi}{\partial x^2} = i \hbar \frac{\partial \psi}{\partial t}$

// let us define the terms in natural units such that: $\hbar = 2m = 1$
// and let us ignore iota in the calculations since it is only representing that the time axis is perpendicular to the spatial axis.

// so we end up with the equation: // $\frac{\partial ^2\psi}{\partial x^2} =  \frac{\partial \psi}{\partial t}$

// using finite difference method, we now write this equation as:
// $\frac{\psi_{i+1}^{x} + \psi_{i-1}^{x} - 2\psi_{i}^{x}}{(\Delta x^2)} = \frac{\psi_{i+1}^{t} - \psi_{i}^{t}}{\Delta t}$

// we will assume that the value of wave function is zero at the boundarys and will take the value of wavefunction at the mid of spatial axis as, say, 5, so as to mimic a spike which will now propagate in spatial domain as the time pases.

n = 50 // no. of points

dx = 0.01 // step size in spatial axis
dt = 0.05 // step size in time axis

hbar = 1.0545e-34  //reduced planck function
m = 9.1e-31 //mass of electron



L = 10 // length of the spatial domain
x = linspace(-L/2, L/2, n) // dividing the x axis in n terms

M = zeros(n,n) // defining a matrix in which as we move along the row, we are moving in spatial axis and as we move along the column, we move in the time axis

M(1,:) = sin(x) + cos(x)
M(:,1) = 0
M(:,n) = 0


for i = 2:(n)
    for j = 2:(n-1)
        M(i,j) = ((dt/((dx)^2)) * (M(i-1,j-1) + M(i-1,j+1) - 2*M(i-1,j))) + M(i,j-1)
    end
end


for i = 1:n
    plot(x,(M(i,:).^2)/1e200)
end
