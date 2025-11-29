using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;

namespace Moria.ReFac
{
    public class GH_SosiImport : GH_Component, IGH_VariableParameterComponent
    {
        public GH_SosiImport() : base(
            "SOSI.Import",
            "SOSI",
            "Leser en SOSI-fil (.sos) og lager sirkler for PUNKT og rør for KURVE.",
            "Tunnel",
            "Import")
        { }

        public override Guid ComponentGuid => new Guid("2E2C3D7E-5C8E-4A19-9D3F-0E6F1C8E7C12");
        protected override System.Drawing.Bitmap Icon => null;

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddTextParameter("SOSI-fil", "path", "Full sti til SOSI-fil (.sos)", GH_ParamAccess.item);
            p.AddNumberParameter("Default punkt-radius [m]", "rP", "Standard radius for PUNKT-sirkler", GH_ParamAccess.item, 1.0);
            p.AddNumberParameter("Default rør-diameter [m]", "dK", "Standard diameter for KURVE-rør", GH_ParamAccess.item, 1.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddCurveParameter("Punkt-sirkler", "Circles", "Sirkler (Curve) i Z=H for alle PUNKT.", GH_ParamAccess.list);
            p.AddPointParameter("Punkt-sentre", "Pts", "Senterpunkt for alle PUNKT.", GH_ParamAccess.list);
            p.AddBrepParameter("Rør", "Pipes", "Rør (Brep) laget fra KURVE-objekter.", GH_ParamAccess.list);
            p.AddCurveParameter("Kurve-midter", "Paths", "Polyline/kurver brukt som senterlinje for rør.", GH_ParamAccess.list);
            p.AddTextParameter("Info", "i", "Status og diagnostikk.", GH_ParamAccess.list);
            p.AddTextParameter("Typer", "types", "Oversikt over funnede OBJTYPE og antall (PUNKT/KURVE).", GH_ParamAccess.list);
        }

        private HashSet<string> _lastPointTypes = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        private HashSet<string> _lastCurveTypes = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

        protected override void SolveInstance(IGH_DataAccess da)
        {
            string path = null;
            double defaultPointRadius = 1.0;
            double defaultPipeDiameter = 1.0;

            if (!da.GetData(0, ref path) || string.IsNullOrWhiteSpace(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler sti til SOSI-fil.");
                return;
            }
            da.GetData(1, ref defaultPointRadius);
            da.GetData(2, ref defaultPipeDiameter);

            if (!File.Exists(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Finner ikke fil: {path}");
                return;
            }

            string[] lines;
            try { lines = File.ReadAllLines(path, Encoding.UTF8); }
            catch (Exception ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Kunne ikke lese fil: {ex.Message}");
                return;
            }

            // --- Headerverdier (fra .HODE-blokka) ---
            double enhet = 1.0;
            double origoN = 0.0;
            double origoO = 0.0;
            bool inHode = false;

            var inv = CultureInfo.InvariantCulture;
            var punkter = new List<SosiPoint>();
            var kurver = new List<SosiCurve>();
            SosiCurve currentCurve = null;
            SosiPoint currentPoint = null;

            foreach (var raw in lines)
            {
                var line = raw.Trim();
                if (line.Length == 0) continue;

                // --- Detekter .HODE-blokka ---
                if (line.StartsWith(".HODE", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = true;
                    continue;
                }

                // Slutt på header når første geometri starter
                if (line.StartsWith(".KURVE", StringComparison.OrdinalIgnoreCase) || line.StartsWith(".PUNKT", StringComparison.OrdinalIgnoreCase))
                    inHode = false;

                // --- LES HEADER INNE I .HODE ---
                if (inHode)
                {
                    if (line.IndexOf("ORIGO", StringComparison.OrdinalIgnoreCase) >= 0)
                    {
                        var nums = SplitNumbers(line);
                        if (nums.Count >= 2)
                        {
                            if (TryParseDoubleInvariant(nums[0], out var n)) origoN = n;
                            if (TryParseDoubleInvariant(nums[1], out var o)) origoO = o;
                        }
                        continue;
                    }

                    if (line.IndexOf("ENHET", StringComparison.OrdinalIgnoreCase) >= 0)
                    {
                        var nums = SplitNumbers(line);
                        if (nums.Count >= 1 && TryParseDoubleInvariant(nums[0], out var scale))
                            enhet = scale;
                        continue;
                    }

                    continue;
                }

                // --- GEOMETRI (etter .HODE) ---
                if (line.StartsWith(".KURVE", StringComparison.OrdinalIgnoreCase))
                {
                    currentCurve = new SosiCurve();
                    currentPoint = null;
                    kurver.Add(currentCurve);
                    continue;
                }
                if (line.StartsWith(".PUNKT", StringComparison.OrdinalIgnoreCase))
                {
                    currentPoint = new SosiPoint();
                    currentCurve = null;
                    punkter.Add(currentPoint);
                    continue;
                }

                if (line.StartsWith("..OBJTYPE", StringComparison.OrdinalIgnoreCase))
                {
                    string val = line.Substring("..OBJTYPE".Length).Trim().Trim(':').Trim();
                    if (currentCurve != null) currentCurve.ObjType = val;
                    if (currentPoint != null) currentPoint.ObjType = val;
                    continue;
                }

                if (line.StartsWith("..NØH", StringComparison.OrdinalIgnoreCase))
                    continue;

                var nums2 = SplitNumbers(line);
                if (nums2.Count >= 3 &&
                    double.TryParse(nums2[0], NumberStyles.Float, inv, out double N) &&
                    double.TryParse(nums2[1], NumberStyles.Float, inv, out double O) &&
                    double.TryParse(nums2[2], NumberStyles.Float, inv, out double H))
                {
                    var pt = new Point3d(
                        origoO + O * enhet,   // X = Ø
                        origoN + N * enhet,   // Y = N
                        H * enhet);           // Z = H

                    if (currentCurve != null) currentCurve.Points.Add(pt);
                    else if (currentPoint != null) currentPoint.Position = pt;
                }
            }

            // --- Typer og dynamiske parametre ---
            var pointTypes = new HashSet<string>(
                punkter.Where(p => p.HasPosition).Select(p => string.IsNullOrWhiteSpace(p.ObjType) ? "(ukjent)" : p.ObjType),
                StringComparer.OrdinalIgnoreCase);

            var curveTypes = new HashSet<string>(
                kurver.Where(k => k.Points.Count >= 2).Select(k => string.IsNullOrWhiteSpace(k.ObjType) ? "(ukjent)" : k.ObjType),
                StringComparer.OrdinalIgnoreCase);

            if (!pointTypes.SetEquals(_lastPointTypes) || !curveTypes.SetEquals(_lastCurveTypes))
            {
                _lastPointTypes = pointTypes;
                _lastCurveTypes = curveTypes;
                EnsureTypeParameters(pointTypes, curveTypes, defaultPointRadius, defaultPipeDiameter);
                ExpireSolution(true);
                return;
            }

            var pointTypeValues = GetDynamicValues(pointTypes, true, defaultPointRadius);
            var curveTypeValues = GetDynamicValues(curveTypes, false, defaultPipeDiameter);

            var infoLines = new List<string>();
            var outCircles = new List<Curve>();
            var outCenters = new List<Point3d>();
            var outPipes = new List<Brep>();
            var outCurves = new List<Curve>();
            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;

            foreach (var pnt in punkter.Where(p => p.HasPosition))
            {
                var t = string.IsNullOrWhiteSpace(pnt.ObjType) ? "(ukjent)" : pnt.ObjType;
                double r = pointTypeValues.TryGetValue(t, out var rv) ? rv : defaultPointRadius;
                if (r <= 0) r = defaultPointRadius;
                var c = new Circle(new Plane(pnt.Position, Vector3d.ZAxis), r);
                outCircles.Add(new ArcCurve(c));
                outCenters.Add(pnt.Position);
            }

            foreach (var k in kurver.Where(k => k.Points.Count >= 2))
            {
                var t = string.IsNullOrWhiteSpace(k.ObjType) ? "(ukjent)" : k.ObjType;
                double d = curveTypeValues.TryGetValue(t, out var dv) ? dv : defaultPipeDiameter;
                if (d <= 0) d = defaultPipeDiameter;
                double r = Math.Max(1e-6, d * 0.5);
                Curve center = new PolylineCurve(k.Points);
                outCurves.Add(center);
                var pipes = Brep.CreatePipe(center, r, false, PipeCapMode.Round, true, tol, RhinoMath.ToRadians(1.0));
                if (pipes != null && pipes.Length > 0)
                    outPipes.Add(pipes[0]);
            }

            infoLines.Add($"Lest: {punkter.Count} PUNKT, {kurver.Count} KURVE.");
            infoLines.Add($"ENHET={enhet}, ORIGO-NØ=({origoN}, {origoO}).");
            infoLines.Add($"DEBUG First KURVE pt (Rhino XYZ): {kurver.FirstOrDefault()?.Points.FirstOrDefault().ToString() ?? "ingen"}");

            da.SetDataList(0, outCircles);
            da.SetDataList(1, outCenters);
            da.SetDataList(2, outPipes);
            da.SetDataList(3, outCurves);
            da.SetDataList(4, infoLines);
            da.SetDataList(5, pointTypes.Concat(curveTypes).ToList());
        }

        // ---------- Hjelpefunksjoner ----------
        private static bool TryParseDoubleInvariant(string s, out double v)
        {
            return double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out v)
                || double.TryParse(s.Replace(',', '.'), NumberStyles.Float, CultureInfo.InvariantCulture, out v);
        }

        // ⚡ Ny enkel SplitNumbers — splitter kun på mellomrom ⚡
        private static List<string> SplitNumbers(string line)
        {
            var res = new List<string>();
            if (string.IsNullOrWhiteSpace(line))
                return res;

            var parts = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
            foreach (var p in parts)
            {
                var s = p.Replace(',', '.');
                if (double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out _))
                    res.Add(s);
            }

            return res;
        }

        // ---------- Dynamiske verdier ----------
        private Dictionary<string, double> GetDynamicValues(HashSet<string> types, bool isPoint, double fallback)
        {
            var map = new Dictionary<string, double>(StringComparer.OrdinalIgnoreCase);
            foreach (var t in types)
            {
                string key = (isPoint ? "P:" : "K:") + t;
                double v = fallback;
                var param = Params.Input.FirstOrDefault(ip => ip.Name.Equals(key, StringComparison.OrdinalIgnoreCase)) as Param_Number;
                if (param != null)
                {
                    if (!TryGetNumber(param, out v))
                        if (!TryGetPersistentDouble(param, out v)) v = fallback;
                }
                map[t] = v;
            }
            return map;
        }

        private bool TryGetPersistentDouble(Param_Number p, out double value)
        {
            value = double.NaN;
            try
            {
                if (p.PersistentDataCount > 0)
                {
                    var branch = p.PersistentData.get_Branch(0);
                    if (branch != null && branch.Count > 0)
                    {
                        var goo = branch[0];
                        if (GH_Convert.ToDouble(goo, out double v, GH_Conversion.Both))
                        {
                            value = v;
                            return true;
                        }
                    }
                }
            }
            catch { }
            return false;
        }

        private bool TryGetNumber(Param_Number p, out double value)
        {
            value = double.NaN;
            if (p.VolatileDataCount > 0)
            {
                var branch = p.VolatileData.get_Branch(0);
                if (branch != null && branch.Count > 0)
                {
                    var goo = branch[0];
                    if (GH_Convert.ToDouble(goo, out double v, GH_Conversion.Both))
                    {
                        value = v;
                        return true;
                    }
                }
            }
            return false;
        }

        // ---------- Variable Parameters ----------
        public bool CanInsertParameter(GH_ParameterSide side, int index) => false;
        public bool CanRemoveParameter(GH_ParameterSide side, int index) => side == GH_ParameterSide.Input && index >= 3;
        public IGH_Param CreateParameter(GH_ParameterSide side, int index) => null;
        public bool DestroyParameter(GH_ParameterSide side, int index) => side == GH_ParameterSide.Input && index >= 3;

        public void VariableParameterMaintenance() { }

        private void EnsureTypeParameters(HashSet<string> pTypes, HashSet<string> kTypes, double defR, double defD) { }

        // ---------- Modeller ----------
        private sealed class SosiPoint
        {
            public string ObjType { get; set; } = string.Empty;
            public Point3d Position { get; set; } = Point3d.Unset;
            public bool HasPosition => Position.IsValid;
        }

        private sealed class SosiCurve
        {
            public string ObjType { get; set; } = string.Empty;
            public List<Point3d> Points { get; } = new List<Point3d>();
        }
    }
}
