#pragma once
#include "generic.h"
#include <fstream>

// What to do about with constructors of abstract classes<

class Config {
	public:		
		Config();
		Config(string fileName);
		~Config();

		string	 instancePath;
		instance loadedInstance;

		double timeout; 

	protected:
		void loadInstance(string fileName);
};
