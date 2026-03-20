using System;
using Rhino.Geometry;

namespace Moria.TunnelGeometry.Components
{
    [Serializable]
    public class PullBoxInstance
    {
        public Guid Id { get; set; } = Guid.NewGuid();

        // Placement
        public double Station { get; set; }
        public double StationOffset { get; set; } = 0.0;
        public bool FollowInterval { get; set; } = true;

        // If true: dims/offsets come from GH inputs (sliders)
        public bool UseDefaults { get; set; } = true;

        // Outer dims (m)
        public double Length { get; set; } = 1.253;  // along tunnel (Z)
        public double Width { get; set; } = 0.580;   // side (X)
        public double Height { get; set; } = 0.500;  // vertical (Y)

        // Openings (m)
        // NOTE: Height is measured from ground up (y=0 -> y=OpenHeight)
        // Long-side opening: width along tunnel (Z)
        public double LongOpenWidth { get; set; } = 0.25;
        public double LongOpenHeight { get; set; } = 0.38;

        // Short-side opening (end face): width along X
        public double ShortOpenWidth { get; set; } = 0.25;
        public double ShortOpenHeight { get; set; } = 0.38;

        // Shell thickness (m) - walls + roof (+bottom in this version)
        public double Thickness { get; set; } = 0.06;

        // Offsets
        public double OffsetZ { get; set; } = 0.0;
        public double OffsetXY { get; set; } = 0.0;

        // Optional cached
        public Plane Frame { get; set; } = Plane.Unset;

        public double EffectiveStation => Station + StationOffset;

        public PullBoxInstance Clone() => (PullBoxInstance)MemberwiseClone();
    }

    [Serializable]
    public class PullBoxPlacementInfo
    {
        public Guid Id { get; set; }
        public double Station { get; set; }
        public Plane Frame { get; set; }

        public double Length { get; set; }
        public double Width { get; set; }
        public double Height { get; set; }

        public double LongOpenWidth { get; set; }
        public double LongOpenHeight { get; set; }

        public double ShortOpenWidth { get; set; }
        public double ShortOpenHeight { get; set; }

        public double Thickness { get; set; }

        public double OffsetZ { get; set; }
        public double OffsetXY { get; set; }

        public bool RightSide { get; set; }
    }
}