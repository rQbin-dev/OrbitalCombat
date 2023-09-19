using Godot;
using Godot.Collections;
using Orbit;
using System;
using System.Linq;

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

	[Export] public Marker3D _marker;

	[Export] public Node _Draw3D;

	
	private Array<MeshInstance3D> _points = new();
	private Array<MeshInstance3D> _lines = new();

	public override void _Ready()
	{
		previousPosition = InitialPosition;
		previousVelocity = InitialVelocity;

		_orbit = OrbitalParameters.From(InitialPosition, InitialVelocity, sizeUnit);
		DrawInitialState();
		DrawOrbit();
	}

	private Vector3 previousVelocity;
	private Vector3 previousPosition;

	private double _currentTime = 0;
	public override void _Process(double delta)
	{
		if (previousPosition != InitialPosition || previousVelocity != InitialVelocity)
		{
			_orbit = OrbitalParameters.From(InitialPosition, InitialVelocity, sizeUnit);
			previousVelocity = InitialVelocity;
			previousPosition = InitialPosition;
			DrawOrbit();
		}


		_currentTime += delta;
	}


	private void DrawInitialState()
	{
		if (_marker == null) return;

		_marker.Position = InitialPosition;
		_marker.Scale = InitialVelocity;
	}


	private void DrawOrbit()
	{
		ClearOrbit();

		Array<Vector3> points = new();
		
		// add points
		for (int i = 0; i < LineSteps; i++)
		{
			double time = (double)i / LineSteps;
			Vector3 point_other = _orbit.calculatePointAtTime(time * Math.PI * 2);

			points.Add(point_other);
		}

		for(int i = 0; i < points.Count - 2; i++)
		{
			 MeshInstance3D line = _Draw3D.Call("line", points[i], points[i + 1]).As<MeshInstance3D>();
			_lines.Add(line);
		}


		double t_peri = _orbit.getTimeInOrbit(0);
		double t_apo = _orbit.getTimeInOrbit(Math.PI);

		Vector3 periapsis = _orbit.calculatePointAt(0);
		Vector3 apoapsis = _orbit.calculatePointAt(Math.PI);

		Vector3 periapsis1 = _orbit.calculatePointAtTime(t_peri);
		Vector3 apoapsis1 = _orbit.calculatePointAtTime(t_apo);

		GD.PrintS(apoapsis, ":", apoapsis1);

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
