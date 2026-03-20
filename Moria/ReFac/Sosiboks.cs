using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Moria.ReFac
{
    public class GH_SosiBoks : GH_Component
    {
        public GH_SosiBoks() : base(
            "SosiBoks",
            "SosiBoks",
            "Leser en SOSI-fil og lager en boks rundt hvert PUNKT. Punktet ligger i sentrum av boksen.",
            "Tunnel",
            "Import")
        { }

        public override Guid ComponentGuid => new Guid("A13F4D8C-6C9B-4C1E-9F44-1B7E3A912345");
        protected override System.Drawing.Bitmap Icon => null;

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddTextParameter("SOSI-fil", "path", "Full sti til SOSI-fil (.sos)", GH_ParamAccess.item);
            p.AddNumberParameter("Høyde", "H", "Boksens høyde (Z)", GH_ParamAccess.item, 2.0);
            p.AddNumberParameter("Bredde", "B", "Boksens bredde (X)", GH_ParamAccess.item, 1.0);
            p.AddNumberParameter("Lengde", "L", "Boksens lengde (Y)", GH_ParamAccess.item, 1.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddBrepParameter("Bokser", "Boxes", "Bokser generert fra PUNKT.", GH_ParamAccess.list);
            p.AddPointParameter("Punkter", "Pts", "Leste PUNKT fra SOSI.", GH_ParamAccess.list);
            p.AddTextParameter("Info", "i", "Status og diagnostikk.", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            string path = null;
            double height = 2.0;
            double width = 1.0;
            double length = 1.0;

            if (!da.GetData(0, ref path) || string.IsNullOrWhiteSpace(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mangler sti til SOSI-fil.");
                return;
            }

            da.GetData(1, ref height);
            da.GetData(2, ref width);
            da.GetData(3, ref length);

            if (!File.Exists(path))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Finner ikke fil: {path}");
                return;
            }

            var lines = File.ReadAllLines(path, Encoding.UTF8);

            double enhet = 1.0;
            double origoN = 0.0;
            double origoO = 0.0;
            bool inHode = false;
            bool inPunkt = false;

            var punkter = new List<Point3d>();
            var inv = CultureInfo.InvariantCulture;

            foreach (var raw in lines)
            {
                var line = raw.Trim();
                if (line.Length == 0) continue;

                if (line.StartsWith(".HODE", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = true;
                    continue;
                }

                if (line.StartsWith(".PUNKT", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = false;
                    inPunkt = true;
                    continue;
                }

                if (line.StartsWith(".KURVE", StringComparison.OrdinalIgnoreCase))
                {
                    inPunkt = false;
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

                if (inPunkt)
                {
                    var nums = SplitNumbers(line);
                    if (nums.Count >= 3 &&
                        double.TryParse(nums[0], NumberStyles.Float, inv, out double N) &&
                        double.TryParse(nums[1], NumberStyles.Float, inv, out double O) &&
                        double.TryParse(nums[2], NumberStyles.Float, inv, out double H))
                    {
                        var pt = new Point3d(
                            origoO + O * enhet,
                            origoN + N * enhet,
                            H * enhet);

                        punkter.Add(pt);
                        inPunkt = false;
                    }
                }
            }

            var bokser = new List<Brep>();

            foreach (var pt in punkter)
            {
                // Lag box centered på punkt
                var plane = new Plane(pt, Vector3d.ZAxis);

                var intervalX = new Interval(-width * 0.5, width * 0.5);
                var intervalY = new Interval(-length * 0.5, length * 0.5);
                var intervalZ = new Interval(-height * 0.5, height * 0.5);

                var box = new Box(plane, intervalX, intervalY, intervalZ);
                bokser.Add(box.ToBrep());
            }

            var info = new List<string>
            {
                $"Lest {punkter.Count} PUNKT.",
                $"Generert {bokser.Count} bokser."
            };

            da.SetDataList(0, bokser);
            da.SetDataList(1, punkter);
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
