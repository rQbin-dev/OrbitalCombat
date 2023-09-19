using Godot;
using System;
namespace Orbit
{

	public struct Kepler
	{
		//from Sebastian Lagues video about orbits
		public static double SolveKepler(double meanAnomaly, double eccentricity, int maxIterations = 100)
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
		static double KeplerEquation(double E, double M, double e)
		{
			// Here the equation has been rearranged. We're trying to find the value for E where this will return 0.
			return M - E + e * Math.Sin(E);
		}
	}



	public struct OrbitalParameters
	{
		public double SemiMajorAxis;
		public double Eccentricity;
		public double Inclination;
		public double RightAscensionOfAscendingNode;
		public double ArgumentOfPeriapsis;
		public double TrueAnomaly;

		public OrbitalParameters(double a, double e, double i, double W, double w, double v)
		{
			SemiMajorAxis = a;
			Eccentricity = e;
			Inclination = i;
			RightAscensionOfAscendingNode = W;
			ArgumentOfPeriapsis = w;
			TrueAnomaly = v;
		}

		
		public Vector3 calculatePointInOrbit(double timeSeconds)
		{
			double meanAnomaly = 2.0 * Math.PI * timeSeconds;

			//solve the eccentric anomaly with the newton-raphson method
			double eccentricAnomaly = Kepler.SolveKepler(meanAnomaly, Eccentricity);

			// P Q are the coordinates on the 2d plane of the orbit
			double P = SemiMajorAxis * (Math.Cos(eccentricAnomaly) - Eccentricity);
			double Q = SemiMajorAxis * Math.Sin(eccentricAnomaly) * Math.Sqrt(1 - Math.Pow(Eccentricity, 2));
            // TODO: Remove
            if (timeSeconds == 0)
			{
				GD.PrintErr(P, " , ", Q);
			}

			double x, y, z;

			//return new((float)P, (float)Q, (float)0);

			// rotate by argument of periapsis
			x = Math.Cos(ArgumentOfPeriapsis) * P - Math.Sin(ArgumentOfPeriapsis) * Q;
			y = Math.Sin(ArgumentOfPeriapsis) * P + Math.Cos(ArgumentOfPeriapsis) * Q;

			// TODO: Remove
			if (timeSeconds == 0)
			{
				GD.PrintErr("after argument of periapsis: ", new Vector3((float)x, (float)y, 0));
			}

			// rotate by inclination
			z = Math.Sin(Inclination) * y;
			y = Math.Cos(Inclination) * y;
            // TODO: Remove
            if (timeSeconds == 0)
			{
				GD.PrintErr("after inclination: ", new Vector3((float)x, (float)y, 0));
			}

			if (!double.IsNaN(RightAscensionOfAscendingNode))
			{
				// rotate by longitude of ascending node
				var xtemp = x;
				x = Math.Cos(RightAscensionOfAscendingNode) * xtemp - Math.Sin(RightAscensionOfAscendingNode) * y;
				y = Math.Sin(RightAscensionOfAscendingNode) * xtemp + Math.Cos(RightAscensionOfAscendingNode) * y;
			}
            // TODO: Remove
            if (timeSeconds == 0)
			{
				GD.PrintS(meanAnomaly, eccentricAnomaly, P, Q);
				GD.PrintErr(new Vector3((float)x, (float)y, (float)z));
			}
			

			return new((float)x, (float)y, (float)z);
		}

		

	

		public OrbitalParameters GetOrbitalParameters(Vector3 position, Vector3 velocity, float mu)
		{

            double semiMajorAxis;
            double eccentricityLength;
            double inclination;
            double ascAscendingNode;
            double argumentPeriapsis;
            double trueAnomaly;


            Vector3 angularMomentum = position.Cross(velocity);
            Vector3 normal = Vector3.Back.Cross(angularMomentum);
            double angularMomentumLength = angularMomentum.Length();
            double semiLatusRectum = (angularMomentumLength * angularMomentumLength) / mu;


            inclination = CalculateInclination(angularMomentum, angularMomentumLength);

			ascAscendingNode = CalculateRightAscensionOfAscendingNode(normal, angularMomentum);
		
			Vector3 eccentricity_vec = CalculateEccentricityVector(angularMomentum, position, velocity, mu);
			eccentricityLength = eccentricity_vec.Length();

			semiMajorAxis = CalculateSemiMajorAxis(eccentricityLength, semiLatusRectum);

			argumentPeriapsis = CalculateArgumentOfPeriapsis(normal, eccentricity_vec);

			trueAnomaly = CalculateTrueAnomaly(mu, semiLatusRectum, position, velocity);


			GD.PrintS(semiMajorAxis, eccentricityLength, "(", eccentricity_vec, ")", inclination * (180 / Math.PI), ascAscendingNode * (180 / Math.PI), argumentPeriapsis * (180 / Math.PI), trueAnomaly * (180 / Math.PI));
			
			return new OrbitalParameters(semiMajorAxis, eccentricityLength, inclination, ascAscendingNode, argumentPeriapsis, trueAnomaly);
		}

		private static double CalculateInclination(Vector3 H, double h)
		{
			return Math.Acos(H.Z / h);
        }

		private static double CalculateRightAscensionOfAscendingNode(Vector3 n, Vector3 h)
		{
			if (n.IsEqualApprox(Vector3.Zero)) return double.NaN;
            return (2.0 * Math.PI + Math.Atan2(h.X, -h.Y)) % (2 * Math.PI);
        }

		private static Vector3 CalculateEccentricityVector(Vector3 h, Vector3 r, Vector3 v, float mu)
		{
			return v.Cross(h) / mu - r.Normalized();
		}

		private static double CalculateSemiMajorAxis(double e_length, double p)
		{
			return (e_length == 1) ? double.PositiveInfinity : p / (1 - (e_length * e_length));
        }

		private static double CalculateArgumentOfPeriapsis(Vector3 normal, Vector3 eccentricity)
		{
			float ne_length = normal.Length() * eccentricity.Length();

			double w;
			//on equatorial orbits assume that:
			if (ne_length == 0) 
				w = Math.Acos(eccentricity.X / eccentricity.Length());
			else 
				w = Math.Acos(normal.Dot(eccentricity) / ne_length);

			return (eccentricity.Z <= 0) ? 2.0 * Math.PI - w : w;
		}


		private static double CalculateTrueAnomaly(double mu, double p, Vector3 position, Vector3 velocity)
		{
			return (2.0 * Math.PI + Math.Atan2(Math.Sqrt(p / mu) * velocity.Dot(position), p - position.Length())) % (2 * Math.PI);
		}
	}

}
