#pragma once

#include <functional>
#include "../quadrature/monte-carlo.h"
#include <PerlinNoise/PerlinNoise.hpp>
#include <cimg-all.h>

using Function4D = std::tuple<std::function<double(double,double,double, double)>,std::function<double(double,double,double,double,double,double, double, double)>>;

double cubicInterpolate (double p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double nCubicInterpolate (int n, double* p, double coordinates[]) {
	if (n == 1) {
		return cubicInterpolate(p, *coordinates);
	}
	else {
		double arr[4];
		int skip = 1 << (n - 1) * 2;
		arr[0] = nCubicInterpolate(n - 1, p, coordinates + 1);
		arr[1] = nCubicInterpolate(n - 1, p + skip, coordinates + 1);
		arr[2] = nCubicInterpolate(n - 1, p + 2*skip, coordinates + 1);
		arr[3] = nCubicInterpolate(n - 1, p + 3*skip, coordinates + 1);
		return cubicInterpolate(arr, *coordinates);
	}
}

Function4D function4d(int& i, int argc, char **argv) {
    std::function<double(double,double,double, double)> func = [] (double x, double y, double z, double w) { return 0.5f; };
	std::function<double(double,double,double,double,double,double, double, double)> groundtruth;
    unsigned long ground_truth_samples = 128000;
    std::size_t seed = std::random_device()();
    bool has_groundtruth = false;
    for (int j = 0; j < argc-1; ++j) {
        if (std::string(argv[j])=="-ground-truth-samples") ground_truth_samples = atol(argv[++j]);
        else if (std::string(argv[j])=="-ground-truth-seed") seed = atol(argv[++j]);
    }
    
    if (std::string(argv[i]) == "lightfield") {
        // cimg_library::CImg<double> input(argv[++i]);
        // cimg_library::CImg<double> grayscale = input.get_channel(0);
        // for (int c=1;c<input.spectrum();++c) grayscale += input.get_channel(c);
        // grayscale /= grayscale.max();
        // func = [grayscale] (double x, double y, double z)
        //      { return grayscale.cubic_atXY(
        //         float(x)*float(grayscale.width()),
        //         float(1.0f-y)*float(grayscale.height()),std::max(0.0f,std::min(float(grayscale.depth()-1),float(z)*float(grayscale.depth()))),0); };

        // Store all the images
        std::vector<cimg_library::CImg<double>> lightField;
        // Load all the images
        #ifdef _WIN32
        const std::string SEPARATOR_LF = "\\";
        #else 
        const std::string SEPARATOR_LF = "/";
        #endif
        const int TOTAL_IMAGES_LF = 81;
        const std::string FILE_NAME_LF = "input_Cam";
        const std::string FILE_PATH_LF = argv[++i];
        lightField.reserve(TOTAL_IMAGES_LF);
        for(int id=0; id<TOTAL_IMAGES_LF; id++) {
            std::ostringstream file;
            file << FILE_PATH_LF << SEPARATOR_LF << FILE_NAME_LF;
            file << std::internal << std::setfill('0') << std::setw(3) << id;
            file << ".png";

            cimg_library::CImg<double> input(file.str().c_str());
            cimg_library::CImg<double> grayscale = input.get_channel(0);
            for (int c=1;c<input.spectrum();++c) grayscale += input.get_channel(c);
            grayscale /= grayscale.max();

            lightField.push_back(grayscale);
        }
        // Return function
        func = [lightField] (double fx, double fy, double fz, double fw)
             {  
                fx *= float(lightField[0].width());
                fy = (1.0-fy) * float(lightField[0].height());
                fz *= 8.0;
                fw *= 8.0;
                
                const int
                x = (int)fx - (fx>=0?0:1),
                y = (int)fy - (fy>=0?0:1),
                z = (int)fz - (fz>=0?0:1),
                w = (int)fw - (fw>=0?0:1);

                const float dx = fx - x, dy = fy - y, dz = fz - z, dw = fw - w;

                double data[256];
                int index = 0;
                for(int xx=x-1; xx<=x+2; xx++) {
                    for(int yy=y-1; yy<=y+2; yy++) {
                        for(int zz=z-1; zz<=z+2; zz++) {
                            for(int ww=w-1; ww<=w+2; ww++) {
                                
                                int ww_ = ww;
                                int zz_ = zz;
                                if (ww < 0) ww_=0;
                                if (ww > 8) ww_=8;
                                if (zz < 0) zz_=0;
                                if (zz > 8) zz_=8;
                                
                                //std::cout << "Loading LF " << xx << "," << yy << "," << zz_ << "," << ww_ << std::endl;
                                double pixel = lightField[zz_*9+ww_].atXY(xx,yy);
                                data[index] = pixel;
                                index++;
                            }
                        }
                    }
                }
                
                double coordinates[4] = {dx, dy, dz, dw};
                //double coordinates[4] = {fx, fy, fz, fw};
                return nCubicInterpolate(4, data, coordinates);
            };
    }

    if (!has_groundtruth) groundtruth = [func,ground_truth_samples,seed] (double a0, double a1, double a2, double a3, double b0, double b1, double b2, double b3) {
        return viltrum::integrator_monte_carlo_uniform(ground_truth_samples,seed).integrate([func] (const std::array<double,4>& x) { return func(std::get<0>(x),std::get<1>(x), std::get<2>(x), std::get<3>(x)); },viltrum::range(a0,a1,a2,a3,b0,b1,b2,b3));
	};

	return std::make_tuple(func,groundtruth);     
}

Function4D function4d(int argc, char **argv) {
	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-function") 
            return function4d(++i,argc,argv);
	}
    
    return Function4D(
        [] (double x, double y, double z, double w) { return 0.5f; },
        [] (double ax, double ay, double az, double aw, double bx, double by, double bz, double bw) { return 0.5f*(bx-ax)*(by-ay)*(bz-az)*(bw-aw);});
}

