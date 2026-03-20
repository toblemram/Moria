using System.Collections.Generic;
using Rhino.Geometry;

namespace Moria.ReFac
{
    internal sealed class SosiData
    {
        public double Enhet { get; set; } = 1.0;
        public double OrigoN { get; set; } = 0.0;
        public double OrigoO { get; set; } = 0.0;

        public List<SosiPoint> Points { get; } = new();
        public List<SosiCurve> Curves { get; } = new();
        public List<SosiFlate> Flater { get; } = new();
    }

    internal sealed class SosiPoint
    {
        public string ObjType { get; set; } = "";
        public Point3d Pos { get; set; } = Point3d.Unset;
    }

    internal sealed class SosiCurve
    {
        public string ObjType { get; set; } = "";
        public List<Point3d> Points { get; } = new();
    }

    internal sealed class SosiFlate
    {
        public string ObjType { get; set; } = "";
        public List<Point3d> Points { get; } = new();
    }
}
