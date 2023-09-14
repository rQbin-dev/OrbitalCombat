using Godot;
using System;





public partial class Orbit : Line2D
{

	public const double GRAVITATIONAL_CONSTANT = 6.6743015;

    private Ellipse _ellipse;

	[Export] public int LineSteps = 100;

    [ExportGroup("Starting Orbit")]
    [Export] public Vector2 StartVelocity;
    [Export] public Vector2 StartPosition;
    [Export] public float CentralBodyMass;


    private bool _updateEllipse = false;



    // Called when the node enters the scene tree for the first time.
    public override void _Ready()
	{
        GetEllipseParametersFrom(StartVelocity, StartPosition, CentralBodyMass);
	}

	// Called every frame. 'delta' is the elapsed time since the previous frame.
	public override void _Process(double delta)
	{
        if(_updateEllipse)
        {
            CalculateOrbitPoints();
            _updateEllipse = false;
        }

	}
    private void CalculateOrbitPoints()
    {
        for(int i = 0; i < LineSteps; i++)
        {
            double time = i / LineSteps;
            Vector2 point = _ellipse.CalculatePointOnOrbit(time);
            AddPoint(point, i);
        }
    }

	public void GetEllipseParametersFrom(Vector2 velocity, Vector2 position, float central_body_mass)
	{
        double standardGravitationalParam = GRAVITATIONAL_CONSTANT * central_body_mass;

		double angularMomentum = position.Cross(velocity);
        double radius = position.Length();

        double semiMajorAxis = (standardGravitationalParam * radius) / (2.0 * standardGravitationalParam - radius * (velocity.X * velocity.X + velocity.Y * velocity.Y));

        double eccentric_x = position.X / radius - (angularMomentum * velocity.Y) / standardGravitationalParam;
        double eccentric_y = position.Y / radius - (angularMomentum * velocity.X) / standardGravitationalParam;

        double semiMinorAxis = semiMajorAxis * Math.Sqrt(1 - (Math.Pow(eccentric_x, 2) + Math.Pow(eccentric_y, 2)));

		double eccentricity = Math.Sqrt(eccentric_x * eccentric_x + eccentric_y * eccentric_y);


        if (_ellipse == null)
        {
            _ellipse = new(semiMajorAxis, semiMinorAxis, eccentricity);
        }
        else
        {
            _ellipse.SemiMajorAxis = semiMajorAxis;
            _ellipse.SemiMinorAxis = semiMinorAxis;
            _ellipse.Eccentricity = eccentricity;
        }

        _updateEllipse = true;
    }
}