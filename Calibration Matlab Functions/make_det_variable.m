function [det] = make_det_variable(DetectorWidth, DetectorHeight, PixelSize, Magnification)

det.DetectorWidth = DetectorWidth;
det.DetectorHeight = DetectorHeight;
det.PixelSize = PixelSize;
det.Mag = Magnification;