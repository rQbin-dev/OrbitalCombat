using Godot;
using System;

public partial class Ellipse
{
    public double SemiMajorAxis;
    public double SemiMinorAxis;
    public double Eccentricity;

    public Ellipse(double semiMajorAxis, double semiMinorAxis, double eccentricity)
    {
        SemiMajorAxis = semiMajorAxis;
        SemiMinorAxis = semiMinorAxis;
        Eccentricity = eccentricity;
    }

    public Ellipse(double periapsis, double apoapsis)
    {
        SemiMajorAxis = (periapsis + apoapsis) / 2;
        double linearEccentricity = SemiMajorAxis - apoapsis;
        Eccentricity = linearEccentricity / apoapsis;
        SemiMinorAxis = SemiMajorAxis * Math.Sqrt(1 - Eccentricity * Eccentricity);
    }

    public Vector2 CalculatePointOnOrbit(double delta_time)
    {
        double meanAnomaly = 2 * Math.PI * delta_time;

        //solve the eccentric anomaly with the newton-raphson method
        double eccentricAnomaly = SolveKepler(meanAnomaly, Eccentricity);

        float pointX = (float)(SemiMajorAxis * Math.Cos(eccentricAnomaly));
        float pointY = (float)(SemiMinorAxis * Math.Sin(eccentricAnomaly));

        return new Vector2(pointX, pointY);
    }

    //from Sebastian Lagues video about orbits
    public double SolveKepler(double meanAnomaly, double eccentricity, int maxIterations = 100)
    {
        const double h = 0.0001; // step size for approximating gradient of the function
        const double acceptableError = 0.00000001;
        double guess = meanAnomaly;

        for (int i = 0; i < maxIterations; i++)
        {
            double y = KeplerEquation(guess, meanAnomaly, eccentricity);
            // Exit early if output of function is very close to zero
            if (Math.Abs(y) < acceptableError)
            {
                break;
            }
            // Update guess to value of x where the slope of the function intersects the x-axis
            double slope = (KeplerEquation(guess + h, meanAnomaly, eccentricity) - y) / h;
            double step = y / slope;
            guess -= step;
        }
        return guess;
    }
    double KeplerEquation(double E, double M, double e)
    {
        // Here the equation has been rearranged. We're trying to find the value for E where this will return 0.
        return M - E + e * Math.Sin(E);
    }
}
