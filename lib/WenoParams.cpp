#include <stdlib.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>

#include <stdexcept>

#include "dog_ini.h"

#include "WenoParams.h"


WENOParams wenoParams;


static std::vector<std::pair<WENOParams::WENOVersion, std::string>>
		wva = {{WENOParams::WENOVersion::FD, "FD"}, 
			   {WENOParams::WENOVersion::JS, "JS"},
			   {WENOParams::WENOVersion::Z, "Z"}};

static std::string wvtos(WENOParams::WENOVersion wv){
	for(int i = 0; i < wva.size(); ++i){
		if(wva[i].first == wv)
			return wva[i].second;
	}
	throw std::logic_error("Can't find string expression for WENOVersion wv.  Expand wva in source file to include new wv values.");	
}

static WENOParams::WENOVersion stowv(std::string s){
	for(int i = 0; i < wva.size(); ++i){
		if(wva[i].second == s)
			return wva[i].first;
	}
	throw std::runtime_error("Unknown WENOVersion string " + s);
}



void WENOParams::init(IniDocument& ini_doc){
	using namespace std;
	weno_version = JS;
	epsilon = 1e-06;
	power_param = 2;
	alpha_scaling = 1.1;

	string section_label = "weno";
	IniDocument::Section& ini_sec = ini_doc[section_label];

	string s_weno_version = ini_sec["weno_version"];
	string s_epsilon = ini_sec["epsilon"];
	string s_power_param = ini_sec["power_param"];
	string s_alpha_scaling = ini_sec["alpha_scaling"];

	if(!s_weno_version.empty())
		weno_version = stowv(s_weno_version);
	if(!s_epsilon.empty())
		epsilon = stod(s_epsilon);
	if(!s_power_param.empty())
		power_param = stod(s_power_param);
	if(!s_alpha_scaling.empty())
		alpha_scaling = stod(s_alpha_scaling);

	if(epsilon <= 0)
		throw std::runtime_error("epsilon <= 0.");
	if(alpha_scaling < 1.0)
		throw std::runtime_error("alpha_scaling < 1.0.");
}

void WENOParams::append_wenohelp(const char* filename){
//TODO Implement (if it makes sense).
}