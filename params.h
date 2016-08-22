#ifndef SPH3_params.h
#define SPH3_params.h

#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream> // ifstream
#include <sstream> // stringstream


namespace params_config {
    
    struct data: std::map <std::string, std::string>
    {
        // Here is a little convenience method...
        bool iskey( const std::string& s ) const
        {
            return count( s ) != 0;
        }
    };
    //---------------------------------------------------------------------------
    // The extraction operator reads configuration::data until EOF.
    // Invalid data is ignored.
    //
    std::istream& operator >> ( std::istream& ins, data& d )
    {
        std::string s, key, value;
        
        // For each (key, value) pair in the file
        while (std::getline( ins, s ))
        {
            std::string::size_type begin = s.find_first_not_of( " \f\t\v" );
            
            // Skip blank lines
            if (begin == std::string::npos) continue;
            
            // Skip commentary
            if (std::string( "#;" ).find( s[ begin ] ) != std::string::npos) continue;
            
            // Extract the key value
            std::string::size_type end = s.find( '=', begin );
            key = s.substr( begin, end - begin );
            
            // (No leading or trailing whitespace allowed)
            key.erase( key.find_last_not_of( " \f\t\v" ) + 1 );
            
            // No blank keys allowed
            if (key.empty()) continue;
            
            // Extract the value (no leading or trailing whitespace allowed)
            begin = s.find_first_not_of( " \f\n\r\t\v", end + 1 );
            end   = s.find_last_not_of(  " \f\n\r\t\v" ) + 1;
            
            value = s.substr( begin, end - begin );
            
            // Insert the properly extracted (key, value) pair into the map
            d[ key ] = value;
        }
        
        return ins;
    }
    
    //---------------------------------------------------------------------------
    // The insertion operator writes all configuration::data to stream.
    //
    std::ostream& operator << ( std::ostream& outs, const data& d )
    {
        data::const_iterator iter;
        for (iter = d.begin(); iter != d.end(); iter++)
            outs << iter->first << " = " << iter->second << std::endl;
        return outs;
    }
    
    
    
}

class Params{
public:
    double init_mass;
    double init_density;
    
    int num_cells_hor;
    int num_cells_ver;
    int num_cells;
    double width, height, dx, dy; //dx is the cell width and dy is the cell height
    double search_radius;
    double kernel_radius;
    double star_mass;
    double damping;
    double pressure;
    double star_rad;
    double t_end, dt;
    int numTimeSteps;
    std::string image_name;
    std::string filename;
    Params(std::string _filename);
    
    
    
};

Params::Params(std::string _filename){
    
    
    params_config::data myconfig;
    params_config::data::const_iterator iter;
    filename = _filename;
    std::string filename_1 = filename + ".ini";
    std::ifstream f(filename_1.c_str());
//f.open(filename_1.c_str());
    
    image_name = filename + ".pgm";
    f >> myconfig;
    f.close();
    
    //Mass
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("mass" == iter->first){
            init_mass  = std::stod (iter->second);
        }
    }
    //Density
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("density" == iter->first){
            init_density  = std::stod (iter->second);
        }
    }
    
    
    //Dt
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("dt" == iter->first){
            dt  = std::stod (iter->second);
        }
    }
    
    //End time
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("time" == iter->first){
            t_end  = std::stod (iter->second);
        }
    }
    
    //Horizontal cells
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("horizontal cells" == iter->first){
            num_cells_hor = std::stoi (iter->second);
        }
    }
    
    //Verticle cells
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("vertical cells" == iter->first){
            num_cells_ver = std::stoi (iter->second);
        }
    }
    
    //Height
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("height" == iter->first){
            height  = std::stod (iter->second);
        }
    }
    
    //Width
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("width" == iter->first){
            width  = std::stod (iter->second);
        }
    }
    
    //Search radius
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("search radius" == iter->first){
            search_radius  = std::stod (iter->second);
        }
    }
    
    //Kernel radius
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("kernel width" == iter->first){
            kernel_radius  = std::stod (iter->second);
        }
    }
    
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("pressure" == iter->first){
            pressure  = std::stod (iter->second);
        }
    }
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("damping" == iter->first){
            damping  = std::stod (iter->second);
        }
    }
    
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("radius" == iter->first){
            star_rad = std::stod (iter->second);
        }
    }
    
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("star mass" == iter->first){
            star_mass = std::stod (iter->second);
        }
    }
    
    numTimeSteps = (int) floor(t_end / dt );
    num_cells = num_cells_hor * num_cells_ver;
    dx = width / (double) num_cells_hor;
    dy = height / (double) num_cells_ver;
}


#endif // PARAMS

