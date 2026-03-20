using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Moria.ReFac
{
    public class GH_SosiSkilt : GH_Component
    {
        public GH_SosiSkilt() : base(
            "SOSI.Skilt",
            "SosiSkilt",
            "Henter PUNKT og KURVE fra SOSI og lager vertikale skilt som følger kurveretning.",
            "Tunnel",
            "Import")
        { }

        public override Guid ComponentGuid => new Guid("5B3C12E7-98AF-4B8E-9E34-6D2FA1E44511");
        protected override System.Drawing.Bitmap Icon => null;

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddTextParameter("SOSI-fil", "path", "Full sti til SOSI-fil (.sos)", GH_ParamAccess.item);
            p.AddNumberParameter("Diameter", "D", "Diameter på skilt", GH_ParamAccess.item, 1.0);
            p.AddNumberParameter("Tykkelse", "T", "Tykkelse på skilt", GH_ParamAccess.item, 0.05);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("Skilt", "Signs", "Vertikale skilt tilpasset kurveretning.", GH_ParamAccess.list);
            p.AddPointParameter("Punkter", "Pts", "Leste PUNKT.", GH_ParamAccess.list);
            p.AddCurveParameter("Senterlinje", "Center", "Første KURVE funnet i SOSI.", GH_ParamAccess.item);
            p.AddTextParameter("Info", "i", "Status og diagnostikk.", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            string path = null;
            double diameter = 1.0;
            double thickness = 0.05;

            if (!da.GetData(0, ref path)) return;
            da.GetData(1, ref diameter);
            da.GetData(2, ref thickness);

            if (!File.Exists(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Finner ikke SOSI-fil.");
                return;
            }

            var lines = File.ReadAllLines(path, Encoding.UTF8);

            double enhet = 1.0;
            double origoN = 0.0;
            double origoO = 0.0;
            bool inHode = false;

            var inv = CultureInfo.InvariantCulture;

            var punkter = new List<Point3d>();
            List<Point3d> curvePts = null;
            bool firstCurveFound = false;

            foreach (var raw in lines)
            {
                var line = raw.Trim();
                if (line.Length == 0) continue;

                if (line.StartsWith(".HODE", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = true;
                    continue;
                }

                if (line.StartsWith(".KURVE", StringComparison.OrdinalIgnoreCase) && !firstCurveFound)
                {
                    inHode = false;
                    curvePts = new List<Point3d>();
                    firstCurveFound = true;
                    continue;
                }

                if (line.StartsWith(".PUNKT", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = false;
                    continue;
                }

                if (inHode)
                {
                    if (line.Contains("ORIGO"))
                    {
                        var nums = SplitNumbers(line);
                        if (nums.Count >= 2)
                        {
                            double.TryParse(nums[0], NumberStyles.Float, inv, out origoN);
                            double.TryParse(nums[1], NumberStyles.Float, inv, out origoO);
                        }
                    }

                    if (line.Contains("ENHET"))
                    {
                        var nums = SplitNumbers(line);
                        if (nums.Count >= 1)
                            double.TryParse(nums[0], NumberStyles.Float, inv, out enhet);
                    }

                    continue;
                }

                var nums2 = SplitNumbers(line);
                if (nums2.Count >= 3 &&
                    double.TryParse(nums2[0], NumberStyles.Float, inv, out double N) &&
                    double.TryParse(nums2[1], NumberStyles.Float, inv, out double O) &&
                    double.TryParse(nums2[2], NumberStyles.Float, inv, out double H))
                {
                    var pt = new Point3d(
                        origoO + O * enhet,
                        origoN + N * enhet,
                        H * enhet);

                    if (curvePts != null && firstCurveFound)
                        curvePts.Add(pt);
                    else
                        punkter.Add(pt);
                }
            }

            if (curvePts == null || curvePts.Count < 2)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Fant ingen KURVE.");
                return;
            }

            var centerCurve = new PolylineCurve(curvePts);

            var outBreps = new List<Brep>();

            foreach (var p in punkter)
            {
                // Finn nærmeste parameter på kurven
                if (!centerCurve.ClosestPoint(p, out double t)) continue;

                var tangent = centerCurve.TangentAt(t);
                tangent.Unitize();

                // Lag lokalt plan
                var xAxis = tangent;
                var zAxis = Vector3d.ZAxis;
                var yAxis = Vector3d.CrossProduct(zAxis, xAxis);
                yAxis.Unitize();

                var plane = new Plane(p, xAxis, zAxis);

                // Sirkel i planet
                var circle = new Circle(plane, diameter * 0.5);
                var curve = new ArcCurve(circle);

                // Ekstruder langs Y (tykkelse)
                var extrusion = Surface.CreateExtrusion(curve, yAxis * thickness);

                if (extrusion != null)
                {
                    var brep = extrusion.ToBrep();
                    if (brep != null)
                        outBreps.Add(brep);
                }
            }

            da.SetDataList(0, outBreps);
            da.SetDataList(1, punkter);
            da.SetData(2, centerCurve);
            da.SetDataList(3, new List<string> { $"Punkt: {punkter.Count}", $"Skilt laget: {outBreps.Count}" });
        }

        private static List<string> SplitNumbers(string line)
        {
            var res = new List<string>();
            var parts = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
            foreach (var p in parts)
            {
                var s = p.Replace(',', '.');
                if (double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out _))
                    res.Add(s);
            }
            return res;
        }
    }
}
