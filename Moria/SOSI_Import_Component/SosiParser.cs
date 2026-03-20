using System;
using System.Collections.Generic;
using System.Globalization;
using Rhino.Geometry;

namespace Moria.ReFac
{
    internal static class SosiParser
    {
        public static SosiData Parse(string[] lines)
        {
            var d = new SosiData();
            var inv = CultureInfo.InvariantCulture;

            bool inHode = false;

            SosiPoint currentPoint = null;
            SosiCurve currentCurve = null;
            SosiFlate currentFlate = null;

            foreach (var raw in lines)
            {
                var line = raw.Trim();
                if (line.Length == 0) continue;

                if (line.StartsWith(".HODE", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = true;
                    currentPoint = null;
                    currentCurve = null;
                    currentFlate = null;
                    continue;
                }

                if (line.StartsWith(".PUNKT", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = false;
                    currentPoint = new SosiPoint();
                    currentCurve = null;
                    currentFlate = null;
                    d.Points.Add(currentPoint);
                    continue;
                }

                if (line.StartsWith(".KURVE", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = false;
                    currentCurve = new SosiCurve();
                    currentPoint = null;
                    currentFlate = null;
                    d.Curves.Add(currentCurve);
                    continue;
                }

                if (line.StartsWith(".FLATE", StringComparison.OrdinalIgnoreCase))
                {
                    inHode = false;
                    currentFlate = new SosiFlate();
                    currentPoint = null;
                    currentCurve = null;
                    d.Flater.Add(currentFlate);
                    continue;
                }

                if (inHode)
                {
                    if (line.IndexOf("ORIGO", StringComparison.OrdinalIgnoreCase) >= 0)
                    {
                        var nums = SplitNumbers(line);
                        if (nums.Count >= 2)
                        {
                            if (TryParseDouble(nums[0], out var n)) d.OrigoN = n;
                            if (TryParseDouble(nums[1], out var o)) d.OrigoO = o;
                        }
                        continue;
                    }

                    if (line.IndexOf("ENHET", StringComparison.OrdinalIgnoreCase) >= 0)
                    {
                        var nums = SplitNumbers(line);
                        if (nums.Count >= 1 && TryParseDouble(nums[0], out var e)) d.Enhet = e;
                        continue;
                    }

                    continue;
                }

                if (line.StartsWith("..OBJTYPE", StringComparison.OrdinalIgnoreCase))
                {
                    var val = line.Substring("..OBJTYPE".Length).Trim().Trim(':').Trim();
                    if (currentPoint != null) currentPoint.ObjType = val;
                    if (currentCurve != null) currentCurve.ObjType = val;
                    if (currentFlate != null) currentFlate.ObjType = val;
                    continue;
                }

                if (line.StartsWith("..NØH", StringComparison.OrdinalIgnoreCase))
                    continue;

                // Koordinater: N O H
                var nums2 = SplitNumbers(line);
                if (nums2.Count >= 3 &&
                    double.TryParse(nums2[0], NumberStyles.Float, inv, out double N) &&
                    double.TryParse(nums2[1], NumberStyles.Float, inv, out double O) &&
                    double.TryParse(nums2[2], NumberStyles.Float, inv, out double H))
                {
                    var pt = new Point3d(
                        d.OrigoO + O * d.Enhet,
                        d.OrigoN + N * d.Enhet,
                        H * d.Enhet);

                    if (currentPoint != null)
                    {
                        if (!currentPoint.Pos.IsValid)
                            currentPoint.Pos = pt;
                    }
                    else if (currentCurve != null)
                    {
                        currentCurve.Points.Add(pt);
                    }
                    else if (currentFlate != null)
                    {
                        currentFlate.Points.Add(pt);
                    }
                }
            }

            d.Points.RemoveAll(p => !p.Pos.IsValid);
            return d;
        }

        private static bool TryParseDouble(string s, out double v)
        {
            return double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out v)
                || double.TryParse(s.Replace(',', '.'), NumberStyles.Float, CultureInfo.InvariantCulture, out v);
        }

        private static List<string> SplitNumbers(string line)
        {
            var res = new List<string>();
            if (string.IsNullOrWhiteSpace(line)) return res;

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
