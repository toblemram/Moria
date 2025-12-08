using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;

using Grasshopper.Kernel;
using Rhino.Geometry;

using Xbim.Common;
using Xbim.IO;
using Xbim.Ifc;

// IFC 2x3 namespaces
using Xbim.Ifc2x3.Kernel;
using Xbim.Ifc2x3.ProductExtension;
using Xbim.Ifc2x3.GeometryResource;
using Xbim.Ifc2x3.TopologyResource;
using Xbim.Ifc2x3.GeometricModelResource;
using Xbim.Ifc2x3.RepresentationResource;
using Xbim.Ifc2x3.PresentationResource;
using Xbim.Ifc2x3.PresentationAppearanceResource;

namespace Moria.TunnelGeometry.Components
{
    public class GH_ExportIfc_Xbim : GH_Component
    {
        public GH_ExportIfc_Xbim()
            : base("Export IFC2x3 (xBIM)", "ExportIFC",
                  "Exports Breps to IFC2x3 using xBIM 5.1 (supports real BREP).",
                  "Tunnel", "Export")
        { }

        // --------------------------------------------------------------------
        protected override void RegisterInputParams(GH_InputParamManager p)
        {
            p.AddBrepParameter("Breps", "B", "Breps to export as IFC2x3 Brep", GH_ParamAccess.list);
            p.AddTextParameter("FileName", "F", "IFC file name", GH_ParamAccess.item, "Tunnel.ifc");
            p.AddTextParameter("Folder", "Dir", "Export folder", GH_ParamAccess.item, "");
            p.AddColourParameter("Color", "C", "Surface color", GH_ParamAccess.item, Color.LightGray);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager p)
        {
            p.AddTextParameter("Info", "I", "Export messages", GH_ParamAccess.list);
        }

        // --------------------------------------------------------------------
        protected override void SolveInstance(IGH_DataAccess da)
        {
            var breps = new List<Brep>();
            string fileName = "";
            string folder = "";
            Color color = Color.LightGray;

            da.GetDataList(0, breps);
            da.GetData(1, ref fileName);
            da.GetData(2, ref folder);
            da.GetData(3, ref color);

            var info = new List<string>();

            if (!Directory.Exists(folder))
            {
                info.Add("Folder does not exist.");
                da.SetDataList(0, info);
                return;
            }

            if (!fileName.EndsWith(".ifc", StringComparison.OrdinalIgnoreCase))
                fileName += ".ifc";

            var path = Path.Combine(folder, fileName);

            // --------------------------------------------------------------------
            // xBIM PROJECT SETUP
            // --------------------------------------------------------------------
            var creds = new XbimEditorCredentials
            {
                ApplicationDevelopersName = "Moria",
                ApplicationFullName = "Moria IFC Exporter",
                ApplicationIdentifier = "MORIA",
                ApplicationVersion = "1.0",
                EditorsOrganisationName = "Moria"
            };

            using (var model = IfcStore.Create(creds, Xbim.Common.Step21.XbimSchemaVersion.Ifc2X3, XbimStoreType.InMemoryModel))
            using (var txn = model.BeginTransaction("IFC Export"))
            {
                // PROJECT
                var project = model.Instances.New<IfcProject>(p =>
                {
                    p.Name = "Moria IFC Export";
                    p.Initialize(ProjectUnits.SIUnitsUK);
                });

                // CONTEXT
                var geomContext = model.Instances
                    .OfType<IfcGeometricRepresentationContext>()
                    .FirstOrDefault();

                // SITE • BUILDING • STOREY
                var site = model.Instances.New<IfcSite>(s => s.Name = "Site");
                var building = model.Instances.New<IfcBuilding>(b => b.Name = "Building");
                var storey = model.Instances.New<IfcBuildingStorey>(s => s.Name = "Storey 1");

                // AGGREGATION RELATIONSHIPS
                model.Instances.New<IfcRelAggregates>(rel =>
                {
                    rel.RelatingObject = project;
                    rel.RelatedObjects.Add(site);
                });

                model.Instances.New<IfcRelAggregates>(rel =>
                {
                    rel.RelatingObject = site;
                    rel.RelatedObjects.Add(building);
                });

                model.Instances.New<IfcRelAggregates>(rel =>
                {
                    rel.RelatingObject = building;
                    rel.RelatedObjects.Add(storey);
                });

                // --------------------------------------------------------------------
                // IFC COLOR STYLE
                // --------------------------------------------------------------------
                var ifcColor = model.Instances.New<IfcColourRgb>(c =>
                {
                    c.Red = color.R / 255.0;
                    c.Green = color.G / 255.0;
                    c.Blue = color.B / 255.0;
                });

                var shading = model.Instances.New<IfcSurfaceStyleShading>(s => s.SurfaceColour = ifcColor);

                var surfaceStyle = model.Instances.New<IfcSurfaceStyle>(ss =>
                {
                    ss.Name = "SurfaceColor";
                    ss.Side = IfcSurfaceSide.BOTH;
                    ss.Styles.Add(shading);
                });

                var styleAssignment = model.Instances.New<IfcPresentationStyleAssignment>(ps =>
                {
                    ps.Styles.Add(surfaceStyle);
                });

                // --------------------------------------------------------------------
                // EXPORT EACH BREP AS IFC BREP (FACETED)
                // --------------------------------------------------------------------
                var meshParams = MeshingParameters.Default;
                int index = 0;

                foreach (var brep in breps)
                {
                    index++;

                    Mesh[] meshes = Mesh.CreateFromBrep(brep, meshParams);
                    if (meshes == null || meshes.Length == 0)
                    {
                        info.Add($"Brep {index}: unable to mesh.");
                        continue;
                    }

                    var mesh = new Mesh();
                    foreach (var m in meshes)
                        mesh.Append(m);

                    mesh.Vertices.CombineIdentical(true, true);
                    mesh.Normals.ComputeNormals();

                    // --------------------------
                    // BUILD IFC FACES
                    // --------------------------
                    var ifcFaces = new List<IfcFace>();

                    foreach (var f in mesh.Faces)
                    {
                        var vids = f.IsTriangle
                            ? new[] { f.A, f.B, f.C }
                            : new[] { f.A, f.B, f.C, f.D };

                        var points = vids.Select(id =>
                        {
                            var v = mesh.Vertices[id];
                            return model.Instances.New<IfcCartesianPoint>(p => p.SetXYZ(v.X, v.Y, v.Z));
                        }).ToList();

                        var poly = model.Instances.New<IfcPolyLoop>(pl =>
                        {
                            foreach (var p in points)
                                pl.Polygon.Add(p);
                        });

                        var bound = model.Instances.New<IfcFaceOuterBound>(b2 =>
                        {
                            b2.Bound = poly;
                            b2.Orientation = true;
                        });

                        var face = model.Instances.New<IfcFace>(ff =>
                        {
                            ff.Bounds.Add(bound);
                        });

                        ifcFaces.Add(face);
                    }

                    // --------------------------
                    // SHELL + FACETED BREP
                    // --------------------------
                    var shell = model.Instances.New<IfcClosedShell>(cs =>
                    {
                        foreach (var face in ifcFaces)
                            cs.CfsFaces.Add(face);
                    });

                    var facetedBrep = model.Instances.New<IfcFacetedBrep>(br =>
                    {
                        br.Outer = shell;
                    });

                    // --------------------------
                    // REPRESENTATION
                    // --------------------------
                    var shapeRep = model.Instances.New<IfcShapeRepresentation>(sr =>
                    {
                        sr.ContextOfItems = geomContext;
                        sr.RepresentationIdentifier = "Body";
                        sr.RepresentationType = "Brep";
                        sr.Items.Add(facetedBrep);
                    });

                    var prodShape = model.Instances.New<IfcProductDefinitionShape>(sh =>
                    {
                        sh.Representations.Add(shapeRep);
                    });

                    // --------------------------
                    // IFC ELEMENT
                    // --------------------------
                    var proxy = model.Instances.New<IfcBuildingElementProxy>(px =>
                    {
                        px.Name = $"Brep_{index}";
                        px.Representation = prodShape;
                        px.ObjectPlacement = storey.ObjectPlacement;
                    });

                    // --------------------------
                    // APPLY STYLE
                    // --------------------------
                    model.Instances.New<IfcStyledItem>(si =>
                    {
                        si.Item = facetedBrep;
                        si.Styles.Add(styleAssignment);
                    });

                    info.Add($"Exported Brep {index} as IfcFacetedBrep");
                }

                txn.Commit();
                model.SaveAs(path);
            }

            info.Add("Saved IFC2x3 file: " + path);
            da.SetDataList(0, info);
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("17FFB85A-011C-44BF-ABDE-524D06A6623C");
    }
}
