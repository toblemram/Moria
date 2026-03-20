using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;

namespace Moria.ReFac
{
    internal static class SosiGeometry
    {
        internal sealed class BrepResult
        {
            public List<Brep> Breps { get; } = new();
            public List<string> Types { get; } = new();
            public List<string> Info { get; } = new();
        }

        internal sealed class CurveResult
        {
            public List<Curve> Curves { get; } = new();
            public List<string> Types { get; } = new();
            public List<string> Info { get; } = new();
        }

        internal sealed class FlateBrepBuild
        {
            public List<Brep> Breps { get; } = new();
            public int FailCount { get; set; } = 0;
            public Dictionary<int, Brep> FlaterByIndex { get; } = new();
        }

        public static List<Point3d> MakePointList(SosiData d)
            => d.Points.Where(p => p.Pos.IsValid).Select(p => p.Pos).ToList();

        public static List<Curve> MakeCurveList(SosiData d)
            => d.Curves.Where(c => c.Points.Count >= 2).Select(c => (Curve)new PolylineCurve(c.Points)).ToList();

        public static FlateBrepBuild MakeFlateBrepList(SosiData d)
        {
            var outp = new FlateBrepBuild();
            for (int i = 0; i < d.Flater.Count; i++)
            {
                var f = d.Flater[i];
                if (f.Points.Count < 3) continue;

                var brep = TryCreatePlanarBrep(f.Points);
                if (brep != null)
                {
                    outp.Breps.Add(brep);
                    outp.FlaterByIndex[i] = brep;
                }
                else outp.FailCount++;
            }
            return outp;
        }

        public static BrepResult BuildPipes(SosiData d, double diameter)
        {
            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            double radius = Math.Max(1e-6, diameter * 0.5);

            var res = new BrepResult();
            int fail = 0;

            foreach (var c in d.Curves.Where(c => c.Points.Count >= 2))
            {
                var path = (Curve)new PolylineCurve(c.Points);
                var pp = Brep.CreatePipe(path, radius, false, PipeCapMode.Round, true, tol, RhinoMath.ToRadians(1.0));
                if (pp != null && pp.Length > 0 && pp[0] != null)
                {
                    res.Breps.Add(pp[0]);
                    res.Types.Add(SosiUtils.NormalizeType(c.ObjType));
                }
                else fail++;
            }

            res.Info.Add("Mode=Rør");
            res.Info.Add($"KURVE={d.Curves.Count(x => x.Points.Count >= 2)}");
            res.Info.Add($"Rør={res.Breps.Count} (feil={fail})");
            res.Info.Add($"Diameter={diameter}");
            res.Info.Add($"ENHET={d.Enhet}, ORIGO(N,O)=({d.OrigoN},{d.OrigoO})");

            return res;
        }

        public static BrepResult BuildBoxes(SosiData d, double height, double width, double length, Curve centerlineOrNull)
        {
            var res = new BrepResult();

            foreach (var p in d.Points.Where(p => p.Pos.IsValid))
            {
                var plane = ComputeOrientedVerticalPlane(p.Pos, centerlineOrNull);

                var ix = new Interval(-width * 0.5, width * 0.5);
                var iy = new Interval(-length * 0.5, length * 0.5);
                var iz = new Interval(-height * 0.5, height * 0.5);

                var box = new Box(plane, ix, iy, iz).ToBrep();
                if (box != null)
                {
                    res.Breps.Add(box);
                    res.Types.Add(SosiUtils.NormalizeType(p.ObjType));
                }
            }

            res.Info.Add("Mode=Boks");
            res.Info.Add($"PUNKT={d.Points.Count(x => x.Pos.IsValid)}");
            res.Info.Add($"Bokser={res.Breps.Count}");
            res.Info.Add($"Dim(H,B,L)=({height},{width},{length})");
            res.Info.Add(centerlineOrNull != null ? "Senterlinje=JA (rotert)" : "Senterlinje=NEI");
            res.Info.Add($"ENHET={d.Enhet}, ORIGO(N,O)=({d.OrigoN},{d.OrigoO})");

            return res;
        }

        public static CurveResult BuildSigns(SosiData d, double radius, Curve centerlineOrNull)
        {
            var res = new CurveResult();
            radius = Math.Max(1e-6, radius);

            int oriented = 0;

            foreach (var p in d.Points.Where(p => p.Pos.IsValid))
            {
                var tan = GetTangentXYAtClosest(p.Pos, centerlineOrNull, out bool usedCenter);
                if (usedCenter) oriented++;

                var plane = MakeVerticalSignPlane(p.Pos, tan);
                var circle = new Circle(plane, radius);
                res.Curves.Add(new ArcCurve(circle));
                res.Types.Add(SosiUtils.NormalizeType(p.ObjType));
            }

            res.Info.Add("Mode=Skilt");
            res.Info.Add($"PUNKT={d.Points.Count(x => x.Pos.IsValid)}");
            res.Info.Add($"Skilt={res.Curves.Count}");
            res.Info.Add($"Radius={radius}");
            res.Info.Add(centerlineOrNull != null ? $"Senterlinje=JA (orientert={oriented})" : "Senterlinje=NEI");
            res.Info.Add($"ENHET={d.Enhet}, ORIGO(N,O)=({d.OrigoN},{d.OrigoO})");

            return res;
        }

        public static BrepResult BuildVifteCylinders(SosiData d, double radius, double dist, Curve centerlineOrNull)
        {
            var res = new BrepResult();
            radius = Math.Max(1e-6, radius);
            dist = Math.Abs(dist);
            if (dist < 1e-6) dist = 0.01;

            int fail = 0;
            int oriented = 0;

            foreach (var p in d.Points.Where(p => p.Pos.IsValid))
            {
                var tan = GetTangentXYAtClosest(p.Pos, centerlineOrNull, out bool usedCenter);
                if (usedCenter) oriented++;

                var circleCenterBase = p.Pos + Vector3d.ZAxis * radius;

                var center0 = circleCenterBase - tan * dist;
                var center1 = circleCenterBase + tan * dist;

                var plane0 = MakeVerticalSignPlane(center0, tan);
                var plane1 = MakeVerticalSignPlane(center1, tan);

                var c0 = new Circle(plane0, radius);
                var c1 = new Circle(plane1, radius);

                var crv0 = c0.ToNurbsCurve();
                var crv1 = c1.ToNurbsCurve();

                var loft = Brep.CreateFromLoft(
                    new List<Curve> { crv0, crv1 },
                    Point3d.Unset,
                    Point3d.Unset,
                    LoftType.Normal,
                    false);

                if (loft != null && loft.Length > 0 && loft[0] != null)
                {
                    res.Breps.Add(loft[0]);
                    res.Types.Add(SosiUtils.NormalizeType(p.ObjType));
                }
                else
                {
                    fail++;
                }
            }

            res.Info.Add("Mode=Vifte");
            res.Info.Add($"PUNKT={d.Points.Count(x => x.Pos.IsValid)}");
            res.Info.Add($"Sylindre={res.Breps.Count} (feil={fail})");
            res.Info.Add($"Radius={radius}");
            res.Info.Add($"Avstand={dist} (total lengde ~ {2.0 * dist})");
            res.Info.Add(centerlineOrNull != null ? $"Senterlinje=JA (orientert={oriented})" : "Senterlinje=NEI");
            res.Info.Add($"ENHET={d.Enhet}, ORIGO(N,O)=({d.OrigoN},{d.OrigoO})");

            return res;
        }

        public static Dictionary<string, List<object>> BuildObjTypeBuckets(
            SosiData d,
            Dictionary<int, Brep> flaterByFlateIndex,
            List<string> typeKeys)
        {
            var map = new Dictionary<string, List<object>>(StringComparer.OrdinalIgnoreCase);
            foreach (var k in typeKeys) map[k] = new List<object>();

            foreach (var p in d.Points.Where(p => p.Pos.IsValid))
            {
                var t = SosiUtils.NormalizeType(p.ObjType);
                if (!map.ContainsKey(t)) map[t] = new List<object>();
                map[t].Add(new GH_Point(p.Pos));
            }

            foreach (var c in d.Curves.Where(c => c.Points.Count >= 2))
            {
                var t = SosiUtils.NormalizeType(c.ObjType);
                if (!map.ContainsKey(t)) map[t] = new List<object>();
                map[t].Add(new GH_Curve(new PolylineCurve(c.Points)));
            }

            for (int i = 0; i < d.Flater.Count; i++)
            {
                if (!flaterByFlateIndex.TryGetValue(i, out var brep)) continue;

                var t = SosiUtils.NormalizeType(d.Flater[i].ObjType);
                if (!map.ContainsKey(t)) map[t] = new List<object>();
                map[t].Add(new GH_Brep(brep));
            }

            return map;
        }

        private static Brep TryCreatePlanarBrep(List<Point3d> pts)
        {
            try
            {
                if (pts == null || pts.Count < 3) return null;

                var pl = new Polyline(pts);
                if (!pl.IsClosed) pl.Add(pl[0]);

                var crv = new PolylineCurve(pl);
                if (!crv.IsClosed) return null;

                var tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
                var breps = Brep.CreatePlanarBreps(crv, tol);
                return (breps != null && breps.Length > 0) ? breps[0] : null;
            }
            catch { return null; }
        }

        private static Plane ComputeOrientedVerticalPlane(Point3d origin, Curve centerline)
        {
            var tan = GetTangentXYAtClosest(origin, centerline, out _);

            var z = Vector3d.ZAxis;
            var y = tan;
            var x = Vector3d.CrossProduct(z, y);
            if (x.IsTiny(1e-9)) x = Vector3d.XAxis;
            else x.Unitize();

            return new Plane(origin, x, y);
        }

        private static Vector3d GetTangentXYAtClosest(Point3d origin, Curve centerline, out bool usedCenterline)
        {
            usedCenterline = false;

            if (centerline != null && centerline.IsValid && centerline.GetLength() > 0)
            {
                if (centerline.ClosestPoint(origin, out double t))
                {
                    var tan = centerline.TangentAt(t);
                    tan.Z = 0.0;
                    if (!tan.IsTiny(1e-9))
                    {
                        tan.Unitize();
                        usedCenterline = true;
                        return tan;
                    }
                }
            }

            var fallback = Vector3d.YAxis;
            fallback.Unitize();
            return fallback;
        }

        private static Plane MakeVerticalSignPlane(Point3d center, Vector3d tanXY)
        {
            var z = Vector3d.ZAxis;

            var tan = tanXY;
            if (tan.IsTiny(1e-9)) tan = Vector3d.YAxis;
            else tan.Unitize();

            var x = Vector3d.CrossProduct(z, tan);
            if (x.IsTiny(1e-9)) x = Vector3d.XAxis;
            else x.Unitize();

            return new Plane(center, x, z);
        }
    }
}
