#pragma once
#include <cimg-all.h>

template<typename T>
class CImgWrapper {
    cimg_library::CImg<T> image;
public:
    class Pixel {
        cimg_library::CImg<T>& image;
        int i, j; 
    public:
        Pixel(cimg_library::CImg<T>& image, int i, int j) : image(image), i(i), j(j) { }
		template<typename A>
        Pixel& operator=(const A& p) {
            for (int c = 0; c<3; ++c) image(i,j,0,c)=p[c]; 
			return (*this);
        }
		template<typename A>
        Pixel& operator+=(const A& p) {
            for (int c = 0; c<3; ++c) image(i,j,0,c)+=p[c]; 
			return (*this);
        }
		T& operator[](int c) { return image(i,j,0,c); }
    };
    class ConstPixel {
        const cimg_library::CImg<T>& image;
        int i, j; 
    public:
        ConstPixel(const cimg_library::CImg<T>& image, int i, int j) : image(image), i(i), j(j) { }
		const T& operator[](int c) const { return image(i,j,0,c); }

    };
    CImgWrapper(int w, int h) : image(w,h,1,3) { }
    void clear() { image.fill(T(0)); }

    Pixel operator()(const std::array<std::size_t,2>& pos) {
        return Pixel(image,std::get<0>(pos),std::get<1>(pos));
    }
	ConstPixel operator()(const std::array<std::size_t,2>& pos) const {
        return ConstPixel(image,std::get<0>(pos),std::get<1>(pos));
    }

    std::array<std::size_t,2> resolution() const {
        return std::array<std::size_t,2>{std::size_t(image.width()),std::size_t(image.height())};
    }
	
	bool load_hdr(const char* name) {
		try {
			image.load_hdr(name); //General load transformed this into a single channel image which is weird.
			return true;
		} catch (const cimg_library::CImgIOException& e) {
			return false;
		}
	}
	
	bool load_hdr(const std::string& name) {
		return load_hdr(name.c_str());
	}

    void save(const char* name) const {
        image.get_max(T(0.01)).save(name);
    }
	
	void save(const std::string& name) const {
		save(name.c_str());
    }
	
	void print() const {
		image.print();
	}
};


