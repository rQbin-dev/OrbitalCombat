using Godot;
using Godot.Collections;
using Orbit;
using System;
using System.Linq;
[Tool]
public partial class OrbitalPath : Node
{

	private Orbit.OrbitalParameters _orbit;

	[ExportCategory("Orbit Visualization")]
	[Export] public int LineSteps = 100;

	[ExportCategory("Orbit Parameters")]
	[Export] public float sizeUnit = 6700;
	[ExportGroup("Initial Parameters")]
	[Export] public Vector3 InitialPosition;
	[Export] public Vector3 InitialVelocity;

	[Export] public Node3D _body;

	[Export] public Node _Draw3D;

	
	private Array<MeshInstance3D> _points = new();
	private Array<MeshInstance3D> _lines = new();


	private double _currentTime = 0;
	public override void _Ready()
	{
		previousPosition = InitialPosition;
		previousVelocity = InitialVelocity;

		_orbit = OrbitalParameters.From(InitialPosition, InitialVelocity, sizeUnit);


		double ea = 2 * Math.Atan(Math.Sqrt((1 - _orbit.Eccentricity) / (1 + _orbit.Eccentricity)) * Math.Tan(_orbit.TrueAnomaly / 2));

		_currentTime = _orbit.getTimeInOrbit(ea);


		DrawInitialState();
		DrawOrbit();
	}

	private Vector3 previousVelocity;
	private Vector3 previousPosition;

	
	public override void _Process(double delta)
	{
		if (previousPosition != InitialPosition || previousVelocity != InitialVelocity)
		{
			_orbit = OrbitalParameters.From(InitialPosition, InitialVelocity, sizeUnit);
			previousVelocity = InitialVelocity;
			previousPosition = InitialPosition;
			DrawOrbit();
		}


		if (Engine.IsEditorHint()) return;

		_currentTime += delta;

		_body.Position = _orbit.calculatePointAtTime(_currentTime);
	}


	private void DrawInitialState()
	{
		if (_body == null) return;

		_body.Position = InitialPosition;
		//_body.Scale = InitialVelocity;
	}


	private void DrawOrbit()
	{
		ClearOrbit();

		Array<Vector3> points = new();
		
		// add points
		for (int i = 0; i < LineSteps; i++)
		{
			double time = (double)i / LineSteps;
			Vector3 point_other = _orbit.calculatePointAt(time * Math.PI * 2);

			points.Add(point_other);
		}

		for(int i = 0; i < points.Count - 2; i++)
		{
			 MeshInstance3D line = _Draw3D.Call("line", points[i], points[i + 1]).As<MeshInstance3D>();
			_lines.Add(line);
		}



		Vector3 periapsis = _orbit.calculatePointAt(0);
		Vector3 apoapsis = _orbit.calculatePointAt(Math.PI);

		//periapsis
		MeshInstance3D periapsisPoint = _Draw3D.Call("point", periapsis, 0.1, Colors.Blue).As<MeshInstance3D>();
		_points.Add(periapsisPoint);

		//apoapsis
		MeshInstance3D apoapsisPoint = _Draw3D.Call("point", apoapsis, 0.1, Colors.Red).As<MeshInstance3D>();
		_points.Add(apoapsisPoint);
	}

	private void ClearOrbit()
	{
		foreach (var p in _points)
		{
			p.QueueFree();
		}
		_points.Clear();
		foreach (var l in _lines)
		{
			l.QueueFree();
		}
		_lines.Clear();
	}

}
