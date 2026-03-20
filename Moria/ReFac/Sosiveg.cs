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
    public class GH_SosiVeg : GH_Component
    {
        public GH_SosiVeg() : base(
            "SOSI.Sosiveg",
            "Sosiveg",
            "Leser SOSI og lager vei fra KURVE med OBJTYPE = Sosiveg.",
            "Tunnel",
            "Import")
        { }

        public override Guid ComponentGuid => new Guid("7F6B92E3-2A7C-4D14-9E5C-AB91C4EFA321");
        protected override System.Drawing.Bitmap Icon => null;

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddTextParameter("SOSI-fil", "path", "Full sti til SOSI-fil (.sos)", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("Veifl ate", "Road", "Veifl ater laget mellom Sosiveg-kurver.", GH_ParamAccess.list);
            p.AddCurveParameter("Kurver", "Curves", "Senterlinjer fra Sosiveg.", GH_ParamAccess.list);
            p.AddTextParameter("Info", "i", "Status og diagnostikk.", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            string path = null;

            if (!da.GetData(0, ref path) || string.IsNullOrWhiteSpace(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler sti til SOSI-fil.");
                return;
            }

            if (!File.Exists(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Finner ikke fil: {path}");
                return;
            }

            string[] lines;
            try
            {
                lines = File.ReadAllLines(path, Encoding.UTF8);
            }
            catch (Exception ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Kunne ikke lese fil: {ex.Message}");
                return;
            }

            double enhet = 1.0;
            double origoN = 0.0;
            double origoO = 0.0;
            bool inHode = false;

            var inv = CultureInfo.InvariantCulture;

            var sosivegKurver = new List<List<Point3d>>();
            List<Point3d> currentCurve = null;
            bool isSosiveg = false;

            foreach (var raw in lines)
            {
                var line = raw.Trim();
                if (line.Length == 0) continue;

                // HEADER
                if (line.StartsWith(".HODE", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = true;
                    continue;
                }

                if (line.StartsWith(".KURVE", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = false;
                    currentCurve = new List<Point3d>();
                    isSosiveg = false;
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

                // OBJTYPE
                if (line.StartsWith("..OBJTYPE", StringComparison.OrdinalIgnoreCase))
                {
                    var type = line.Substring("..OBJTYPE".Length).Trim().Trim(':').Trim();

                    if (type.Equals("Sosiveg", StringComparison.OrdinalIgnoreCase))
                        isSosiveg = true;

                    continue;
                }

                // Koordinater
                var nums2 = SplitNumbers(line);
                if (nums2.Count >= 3 &&
                    double.TryParse(nums2[0], NumberStyles.Float, inv, out double N) &&
                    double.TryParse(nums2[1], NumberStyles.Float, inv, out double O) &&
                    double.TryParse(nums2[2], NumberStyles.Float, inv, out double H))
                {
                    if (currentCurve == null) continue;

                    var pt = new Point3d(
                        origoO + O * enhet,
                        origoN + N * enhet,
                        H * enhet);

                    currentCurve.Add(pt);
                }

                // Når vi møter ny .KURVE eller slutt, lagre hvis Sosiveg
                if (line.StartsWith(".KURVE", StringComparison.OrdinalIgnoreCase))
                {
                    if (isSosiveg && currentCurve != null && currentCurve.Count >= 2)
                        sosivegKurver.Add(currentCurve);
                }
            }

            // Siste kurve
            if (isSosiveg && currentCurve != null && currentCurve.Count >= 2)
                sosivegKurver.Add(currentCurve);

            var rhinoCurves = sosivegKurver
                .Select(c => new PolylineCurve(c))
                .ToList();

            var outBreps = new List<Brep>();
            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 0.01;

            // Loft mellom nabokurver
            for (int i = 0; i < rhinoCurves.Count - 1; i++)
            {
                var loft = Brep.CreateFromLoft(
                    new List<Curve> { rhinoCurves[i], rhinoCurves[i + 1] },
                    Point3d.Unset,
                    Point3d.Unset,
                    LoftType.Normal,
                    false);

                if (loft != null && loft.Length > 0)
                    outBreps.Add(loft[0]);
            }

            var info = new List<string>
            {
                $"Fant {rhinoCurves.Count} Sosiveg-kurver.",
                $"Genererte {outBreps.Count} veifl ater."
            };

            da.SetDataList(0, outBreps);
            da.SetDataList(1, rhinoCurves);
            da.SetDataList(2, info);
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
