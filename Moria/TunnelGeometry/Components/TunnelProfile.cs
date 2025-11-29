using System;
using System.Collections.Generic;
using Grasshopper.Kernel;

using Rhino;
using Rhino.Geometry;

namespace Moria.TunnelGeometry.Components
{
    public class GH_TunnelProfile : GH_Component
    {
        public GH_TunnelProfile() :
            base("TunnelProfile", "TunnelProfile",
                 "Builds a tunnel cross-section using Norwegian T-profile definitions.",
                 "Tunnel", "Sections")
        { }

        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddTextParameter(
                "ProfileType", "T",
                "Profile type, e.g. T5.5, T7.5, T8.5, T9.5–T14.",
                GH_ParamAccess.item, "T14");

            p.AddCurveParameter(
                "Path", "P",
                "Optional sweep path.", GH_ParamAccess.item);

            p.AddBooleanParameter(
                "LeftToRight", "L2R",
                "Force roof arc orientation left → right.",
                GH_ParamAccess.item, true);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddCurveParameter(
                "Profile", "C",
                "Closed profile polycurve.", GH_ParamAccess.item);

            p.AddCurveParameter(
                "Segments", "S",
                "Profile segments: left Rv, roof Rh, right Rv, bottom.",
                GH_ParamAccess.list);

            p.AddBrepParameter(
                "Tunnel", "B",
                "Swept tunnel volume.", GH_ParamAccess.item);

            p.AddTextParameter(
                "Info", "i",
                "Info / debug text.", GH_ParamAccess.list);

            p.AddGeometryParameter(
                "DebugGeom", "D",
                "Debug circles, points and helper lines in WorldXY.",
                GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess da)
        {
            string type = "T14";
            Curve path = null;
            bool leftToRight = true;

            da.GetData(0, ref type);
            da.GetData(1, ref path);
            da.GetData(2, ref leftToRight);

            type = type.Replace(",", ".").ToUpperInvariant();

            double tol = RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-6;
            var info = new List<string>();
            var debugGeom = new List<GeometryBase>();

            // ---------------- Profile data ----------------
            if (!ProfileType.Profiles.TryGetValue(type, out ProfileType.ProfileParameters par))
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Error,
                    $"Invalid profile type '{type}'. Valid: {string.Join(", ", ProfileType.Profiles.Keys)}");
                return;
            }

            // ---------------- Build profile (WorldXY) ----------------
            PolyCurve profile;
            List<Curve> segments;
            string error;

            bool ok;
            if (ProfileType.IsLowRoof(type))
            {
                ok = LowRoofProfileBuilder.Build(
                    type, par, leftToRight, tol,
                    out profile, out segments, out debugGeom, out error);
            }
            else
            {
                ok = StandardProfileBuilder.Build(
                    par, leftToRight, tol,
                    out profile, out segments, out debugGeom, out error);
            }

            if (!ok)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, error);
                return;
            }

            // ---------------- Orient to path (profile only) ----------------
            if (path != null)
            {
                if (path.PerpendicularFrameAt(path.Domain.T0, out Plane frame))
                {
                    // move midpoint of bottom segment to origin
                    Curve bottom = profile.SegmentCurve(profile.SegmentCount - 1);
                    Point3d b0 = bottom.PointAtStart;
                    Point3d b1 = bottom.PointAtEnd;
                    Point3d midBottom = 0.5 * (b0 + b1);

                    Transform toOrigin = Transform.Translation(
                        -midBottom.X, -midBottom.Y, -midBottom.Z);
                    profile.Transform(toOrigin);

                    Transform orient = Transform.PlaneToPlane(Plane.WorldXY, frame);
                    profile.Transform(orient);
                }
                else
                {
                    AddRuntimeMessage(
                        GH_RuntimeMessageLevel.Warning,
                        "Could not compute PerpendicularFrameAt for Path – profile kept in WorldXY.");
                }
            }

            // ---------------- Sweep ----------------
            Brep swept = SweepAlongPath(profile, path, tol);

            // ---------------- Outputs ----------------
            da.SetData(0, profile);
            da.SetDataList(1, segments);
            da.SetData(2, swept);

            info.Add($"Profile: {type}");
            info.Add($"Yv={par.Yv:0.###}, Rv={par.Rv:0.###}, X={par.X:0.###}, Rh={par.Rh:0.###}");
            info.Add($"Closed={profile.IsClosed}, Sweep={(swept != null)}");
            da.SetDataList(3, info);

            da.SetDataList(4, debugGeom);
        }

        private Brep SweepAlongPath(PolyCurve profile, Curve path, double tol)
        {
            if (path == null)
                return null;

            var sweep = new SweepOneRail
            {
                AngleToleranceRadians = RhinoMath.ToRadians(1.0),
                SweepTolerance = tol
            };

            Brep[] breps = sweep.PerformSweep(path, profile);
            if (breps != null && breps.Length > 0)
                return breps[0];

            AddRuntimeMessage(
                GH_RuntimeMessageLevel.Warning,
                "Sweep failed – check that Path is open/valid and not too kinky.");
            return null;
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                var assembly = System.Reflection.Assembly.GetExecutingAssembly();
                using (var stream = assembly.GetManifestResourceStream("Moria.Resources.TunnelProfile.png"))
                {
                    return new System.Drawing.Bitmap(stream);
                }
            }
        }

        public override Guid ComponentGuid =>
            new Guid("4073A8B0-1C97-46A2-AEA7-4E9DE5CC1B63");
    }
}
