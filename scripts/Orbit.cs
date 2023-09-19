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
		public const float GRAVITANIONAL_CONSTANT = 1.0f;
		/**/
		public double SemiMajorAxis;
		public double Eccentricity;
		public double Inclination;
		public double RightAscensionOfAscendingNode;
		public double ArgumentOfPeriapsis;
		public double TrueAnomaly;
		public double LongitudeOfPeriapsis;

		public OrbitalParameters(double a, double e, double i, double W, double w, double v, double pi)
		{
			SemiMajorAxis = a;
			Eccentricity = e;
			Inclination = i;
			RightAscensionOfAscendingNode = W;
			ArgumentOfPeriapsis = w;
			TrueAnomaly = v;
			LongitudeOfPeriapsis = pi;
		}

		
		public Vector3 calculatePointInOrbit(double timeSeconds)
		{
			double meanAnomaly = 2.0 * Math.PI * timeSeconds;

			//solve the eccentric anomaly with the newton-raphson method
			double eccentricAnomaly = Kepler.SolveKepler(meanAnomaly, Eccentricity);

			// P Q are the coordinates on the 2d plane of the orbit
			double P = SemiMajorAxis * (Math.Cos(eccentricAnomaly) - Eccentricity);
			double Q = SemiMajorAxis * Math.Sin(eccentricAnomaly) * Math.Sqrt(1 - Math.Pow(Eccentricity, 2));

			double x, y, z;

			//return new((float)P, (float)Q, (float)0);

			// rotate by argument of periapsis
			x = Math.Cos(ArgumentOfPeriapsis) * P - Math.Sin(ArgumentOfPeriapsis) * Q;
			y = Math.Sin(ArgumentOfPeriapsis) * P + Math.Cos(ArgumentOfPeriapsis) * Q;
			if (timeSeconds == 0)
			{
				GD.PrintErr("after argument of periapsis: ", new Vector3((float)x, (float)y, 0));
			}

			// rotate by inclination
			z = Math.Sin(Inclination) * y;
			y = Math.Cos(Inclination) * y;

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

			if(timeSeconds == 0)
			{
				GD.PrintS(meanAnomaly, eccentricAnomaly, P, Q);
				GD.PrintErr(new Vector3((float)x, (float)y, (float)z));
			}
			

			return new((float)x, (float)y, (float)z);
		}

		

		public static Vector3 GetPointOnOrbit(double timeSeconds, OrbitalParameters parameters)
		{
			double meanAnomaly = 2.0 * Math.PI * timeSeconds;

			//solve the eccentric anomaly with the newton-raphson method
			double eccentricAnomaly = Kepler.SolveKepler(meanAnomaly, parameters.Eccentricity);

			double nu = 2 * Math.Atan(Math.Sqrt((1 + parameters.Eccentricity) / 1 - parameters.Eccentricity) * Math.Tan(eccentricAnomaly / 2));

			//4
			double r = parameters.SemiMajorAxis * (1 - parameters.Eccentricity * Math.Cos(parameters.TrueAnomaly));
			//5
			double h = Math.Sqrt(GRAVITANIONAL_CONSTANT * parameters.SemiMajorAxis * (1 - (parameters.Eccentricity * parameters.Eccentricity)));
			//6
			double Om = parameters.RightAscensionOfAscendingNode;
			double w = parameters.ArgumentOfPeriapsis;

			double X = r * (Math.Cos(Om) * Math.Cos(w + nu) - Math.Sin(Om) * Math.Sin(w + nu) * Math.Cos(parameters.Inclination));
			double Y = r * (Math.Sin(Om) * Math.Cos(w + nu) + Math.Cos(Om) * Math.Sin(w + nu) * Math.Cos(parameters.Inclination));
			double Z = r * (Math.Sin(parameters.Inclination) * Math.Sin(w + nu));

			return new((float)X, (float)Y, (float)Z);
		}

		public OrbitalParameters GetOrbitalParameters(Vector3 position, Vector3 velocity, float mu)
		{
			double semiMajorAxis, eccentricity, inclination, ascAscendingNode, argumentPeriapsis, trueAnomaly, longitudeOfPeriapsis;

			double angularMomentum;

			double radialVelocity = velocity.Dot(position.Normalized());


			Vector3 angularMomentum_vec = position.Cross(velocity);
			angularMomentum = angularMomentum_vec.Length();


			inclination = Math.Acos(angularMomentum_vec.Z / angularMomentum);

			Vector3 normal = Vector3.Back.Cross(angularMomentum_vec);
			ascAscendingNode = CalculateRightAscensionOfAscendingNode(angularMomentum_vec);

			Vector3 eccentricity_vec = velocity.Cross(angularMomentum_vec) / mu - position / position.Length();
			eccentricity = eccentricity_vec.Length();

			double p = (angularMomentum * angularMomentum) / mu;
			semiMajorAxis = (eccentricity == 1) ? double.PositiveInfinity : p / (1 - (eccentricity * eccentricity));


			argumentPeriapsis = CalculateArgumentOfPeriapsis(normal, eccentricity_vec);


			trueAnomaly = CalculateTrueAnomaly(mu, p, position, velocity);

			longitudeOfPeriapsis = CalculateLongitudeOfPeriapsis(eccentricity_vec, eccentricity, inclination);

			GD.PrintS(semiMajorAxis, eccentricity, inclination * (180 / Math.PI), ascAscendingNode * (180 / Math.PI), argumentPeriapsis * (180 / Math.PI), trueAnomaly * (180 / Math.PI), longitudeOfPeriapsis * (180 / Math.PI));

			return new OrbitalParameters(semiMajorAxis, eccentricity, inclination, ascAscendingNode, argumentPeriapsis, trueAnomaly, longitudeOfPeriapsis);
		}

		private static Vector3 CalculateEccentricityVector(Vector3 position, Vector3 velocity, float mu)
		{
			Vector3 part1 = (velocity.LengthSquared() - mu / position.Length()) * position;
			Vector3 part2 = position.Dot(velocity) * velocity;

			return (part1 - part2) / mu;
		}

		private static double CalculateRightAscensionOfAscendingNode(Vector3 h)
		{
			return (2 * Math.PI + Math.Atan2(h.X, -h.Y)) % (2 * Math.PI);
		}

		private static double CalculateArgumentOfPeriapsis(Vector3 normal, Vector3 eccentricity)
		{
			float ne_length = normal.Length() * eccentricity.Length();

			GD.PrintS(normal, eccentricity, "=> acos", normal.Dot(eccentricity), "/", ne_length, " = ", Math.Acos(normal.Dot(eccentricity) / ne_length) * (180 / Math.PI));
			double w;
			//on equatorial orbits assume that:
			if (ne_length == 0) 
				w = Math.Acos(eccentricity.X / eccentricity.Length());
			else 
				w = Math.Acos(normal.Dot(eccentricity) / ne_length);

			return (eccentricity.Z < 0) ? 2.0 * Math.PI - w : w;
		}


		private static double CalculateTrueAnomaly(double mu, double p, Vector3 position, Vector3 velocity)
		{
			return (2 * Math.PI + Math.Atan2(Math.Sqrt(p / mu) * velocity.Dot(position), p - position.Length())) % (2 * Math.PI);
		}
		private static double CalculateLongitudeOfPeriapsis(Vector3 e, double e_length, double i)
		{

			if (e_length == 0) { return double.NaN; }

			if ((e.Y < 0 && i >= 0 && i <= 90) || (e.Y > 0 && i > 90))
			{
				return 2 * Math.PI - Math.Acos(e.X / e.Length());
			}

			return Math.Acos(e.X / e.Length());
		}
	}

}
