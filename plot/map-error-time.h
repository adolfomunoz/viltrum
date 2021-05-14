#pragma once

#include <vector>
#include <cmath>
#include <svg-cpp-plot/svg-cpp-plot.h>

std::string label(float nsamples) {
	std::stringstream s;
	if (nsamples<0.55)
		s<<"1/"<<std::fixed<<std::setprecision(0)<<1.0f/nsamples;
	else if (nsamples >= 1.e6)
		s<<std::fixed<<std::setprecision(0)<<(nsamples*1.e-6)<<"M";
	else if (nsamples >= 1.e3)
		s<<std::fixed<<std::setprecision(0)<<(nsamples*1.e-3)<<"k";
	else
		s<<std::fixed<<std::setprecision(0)<<nsamples;
	return s.str();
}

svg_cpp_plot::SVG map_error_time(
	const std::vector<int>& mc_spp,
	const std::vector<int>& qdt_its,
	const std::vector<std::vector<float>>& err_,
	const std::vector<std::vector<float>>& time_,
	std::size_t bins, std::string_view cmap, std::size_t DIMBINS = 2) {
		
	std::vector<std::vector<float>> err = err_;
	std::vector<std::vector<float>> time = time_;
	if ((mc_spp[0]==0) && (qdt_its[0]==0)) {
		err[0][0] = 
			std::max(err[1][1],std::max(err[1][0],err[0][1]));
		time[0][0] =
			std::min(time[1][1],std::min(time[1][0],time[0][1]));
	}
	
    float dim_add    = std::pow(3,DIMBINS);
	float dim_factor = std::pow(3,DIMBINS-1)*2;
    
	float width = 200; float spacing = 10;
	std::vector<float> ticks_cvspp;
	std::vector<std::string> labels_cvspp;
	std::size_t di=1;
	while ((qdt_its.size()/di)>6) ++di;
	ticks_cvspp.push_back(0);
	labels_cvspp.push_back("<tspan dy=\"1em\"> </tspan><tspan dy=\"1em\">MC</tspan>");
	for (std::size_t i = 1; i<(qdt_its.size()-di); i+=di) {
		ticks_cvspp.push_back(float(i));
		labels_cvspp.push_back(label(float(dim_add+qdt_its[i]*dim_factor)/bins));
	}
	ticks_cvspp.push_back(float(qdt_its.size()-1));
	labels_cvspp.push_back(label(float(dim_add+qdt_its.back()*dim_factor)/bins));
	
	std::vector<float> ticks_mcspp;
	std::vector<std::string> labels_mcspp;
	di=1;
	while ((mc_spp.size()/di)>12) ++di;
	ticks_mcspp.push_back(0);
	labels_mcspp.push_back("ST");
	for (std::size_t i = 1; i<(mc_spp.size()-di); i+=di) {
		ticks_mcspp.push_back(float(i));
		labels_mcspp.push_back(label(mc_spp[i]));
	}
	ticks_mcspp.push_back(float(mc_spp.size()-1));
	labels_mcspp.push_back(label(mc_spp.back()));

	std::vector<std::vector<float>> time_error(mc_spp.size(),std::vector<float>(qdt_its.size()));
	std::vector<std::vector<float>> samples_error(mc_spp.size(),std::vector<float>(qdt_its.size()));

	for (std::size_t m = 0; m<mc_spp.size(); ++m) {
		for (std::size_t q = 0; q<qdt_its.size(); ++q) {
			time_error[m][q] = std::log(std::exp(time[m][q])*std::exp(err[m][q]));
			samples_error[m][q] =
				std::log( (bins*mc_spp[m]+dim_add + dim_factor*qdt_its[q])*std::exp(err[m][q]) );
		}
	}

	svg_cpp_plot::SVG svg;
   	svg.viewBox(svg_cpp_plot::BoundingBox(
		-40,-40,4*(width+spacing),width+40));
		
	svg_cpp_plot::SVGPlot plt_error;
	plt_error.figsize({width,width}).ylabel("MC spp").xlabel("CV spp").title("Error").xticks(ticks_cvspp).xticklabels(labels_cvspp).yticks(ticks_mcspp).yticklabels(labels_mcspp);
	plt_error.imshow(err).cmap(cmap);
	svg.add(svg_cpp_plot::_2d::group(
		svg_cpp_plot::_2d::translate(0,0))).add(plt_error.graph());
		
	svg_cpp_plot::SVGPlot plt_time;
	plt_time.figsize({width,width}).xlabel("CV spp").title("Time").xticks(ticks_cvspp).xticklabels(labels_cvspp).yticks({});
	plt_time.imshow(time).cmap(cmap);
	svg.add(svg_cpp_plot::_2d::group(
		svg_cpp_plot::_2d::translate(width+spacing,0))).add(plt_time.graph());
	
	svg_cpp_plot::SVGPlot plt_all;
	plt_all.figsize({width,width}).xlabel("CV spp").title("Time X Error").xticks(ticks_cvspp).xticklabels(labels_cvspp).yticks({});
	plt_all.imshow(time_error).cmap(cmap);
	svg.add(svg_cpp_plot::_2d::group(
		svg_cpp_plot::_2d::translate(2*(width+spacing),0))).add(plt_all.graph());
		
	svg_cpp_plot::SVGPlot plt_samples;
	plt_samples.figsize({width,width}).xlabel("CV spp").title("Samples X Error").xticks(ticks_cvspp).xticklabels(labels_cvspp).yticks({});
	plt_samples.imshow(samples_error).cmap(cmap);
	svg.add(svg_cpp_plot::_2d::group(
		svg_cpp_plot::_2d::translate(3*(width+spacing),0))).add(plt_samples.graph());

	return svg;
}