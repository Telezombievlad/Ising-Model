// No Copyright. Vladislav Aleinik 2019

#include "Model.hpp"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <unistd.h>
#include <time.h>
#include <fcntl.h>

#include <sys/mman.h>
#include <sys/ioctl.h>
#include <linux/fb.h>
#include <linux/kd.h>

#include "vendor/cnpy/cnpy.h"

long double getInteractivity(int x, int y)
{
	static std::random_device rd;
	static std::mt19937 gen{rd()};
	static std::uniform_real_distribution<double> distribution(-1.0, 1.0);

	return 1.0;
}

int getStateX(int x, int y)
{
	return rand() % STATE_GRAPH_SIZE_X;
}

int getStateY(int x, int y)
{
	return rand() % STATE_GRAPH_SIZE_Y;
}

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		fprintf(stderr, "[ISING-MODEL] Expected input: model <temperature> <field> <filename>\n");
		exit(EXIT_FAILURE);
	}

	//=================
	// Parse arguments 
	//=================

	double temperature = strtof(argv[1], NULL);
	if (temperature <= 0.0f)
	{
		fprintf(stderr, "[ISING-MODEL] Invalid temperature value: %f\n", temperature);
		exit(EXIT_FAILURE);
	}

	double fieldZ = strtof(argv[2], NULL);

#ifdef RENDERING
	//=================================
	// Open frame buffer for rendering 
	//=================================

	int fb0_fd = open("/dev/fb0", O_RDWR);
	if (fb0_fd == -1)
	{
		fprintf(stderr, "[ISING-MODEL] Unable to open /dev/fb0\n");
		exit(EXIT_FAILURE);
	}
	
	struct fb_var_screeninfo vinf;
    if (ioctl(fb0_fd, FBIOGET_VSCREENINFO, &vinf) == -1)
    {
    	fprintf(stderr, "[ISING-MODEL] Unable get variable screen info\n");
		exit(EXIT_FAILURE);
    }
	
	struct fb_fix_screeninfo finf;
	if (ioctl(fb0_fd, FBIOGET_FSCREENINFO, &finf) == -1)
	{
		fprintf(stderr, "[ISING-MODEL] Unable to get fixed screen info\n");
		exit(EXIT_FAILURE);
	}

	//=========================================
	// Map frame buffer into our address space 
	//=========================================

	char* frame_buffer = (char*) mmap(NULL, finf.line_length * vinf.yres, PROT_READ | PROT_WRITE, MAP_SHARED, fb0_fd, 0);
	if (frame_buffer == MAP_FAILED)
	{
		fprintf(stderr, "[ISING-MODEL] Unable to map frame buffer into address space\n");
		exit(EXIT_FAILURE);
	}
	
	//===================================================
	// Parse frame buffer info (used only for rendering)
	//===================================================

	size_t        fb_size         = finf.line_length * vinf.yres;
	size_t        bytes_per_pixel = vinf.bits_per_pixel/8;
	size_t        bytes_per_line  = finf.line_length;
	uint_fast16_t offset_red      = vinf.red.offset   / 8;
	uint_fast16_t offset_green    = vinf.green.offset / 8;
	uint_fast16_t offset_blue     = vinf.blue.offset  / 8;
	size_t        sizeX           = vinf.xres;
	size_t        sizeY           = vinf.yres;
#else
	size_t sizeX = 2000;
	size_t sizeY = 2000;
#endif 
	
	//===================================
	// Allocate buffer for data triplets 
	//===================================

	size_t saved_data_size = 100;
	size_t save_data_every = 10;
	double* data_triplets = (double*) calloc(3 * saved_data_size, sizeof(*data_triplets));
	if (data_triplets == NULL)
	{
		fprintf(stderr, "[ISING-MODEL] Unable to allocate array for saved data!\n");
		exit(EXIT_FAILURE);
	}

	//==============
	// Calculations 
	//==============

	Lattice isingModel = Lattice(sizeX/2, sizeY/2, getInteractivity, getStateX, getStateY);

	Vector field{0.0, 0.0, fieldZ};
	for (size_t iteration = 0, cur_saved_data = 0;
		 iteration < 100 + save_data_every * saved_data_size;
		 ++iteration)
	{
		for (int i = 0; i < 100000; ++i)
		{
			isingModel.metropolisStep(temperature, field);
		}

		double magnetization = 0.0;
		for (uint_fast16_t x = 0; x < sizeX/2; ++x)
		{
			for (uint_fast16_t y = 0; y < sizeY/2; ++y)
			{
				LatticePoint curPoint = isingModel.get(x, y);
				Vector spin = isingModel.stateGraph.get(curPoint.stateX, curPoint.stateY);

				magnetization += spin.z;

			#ifdef RENDERING
				for (size_t dx = 0; dx < 2; ++dx) {
				for (size_t dy = 0; dy < 2; ++dy) {
					uint_fast16_t pix_offset = (2*x+dx) * bytes_per_pixel +
					                           (2*y+dy) * bytes_per_line;
					frame_buffer[pix_offset + offset_red  ] = 0;
					frame_buffer[pix_offset + offset_green] = (1.0 + spin.z) * 127;
					frame_buffer[pix_offset + offset_blue ] = 32;
				}}
			#endif
			}
		}
		
		magnetization /= sizeX*sizeY/4;

		if (100 <= iteration && iteration % save_data_every == 0 && cur_saved_data < saved_data_size)
		{
			data_triplets[3 * cur_saved_data + 0] = temperature;
			data_triplets[3 * cur_saved_data + 1] = fieldZ;
			data_triplets[3 * cur_saved_data + 2] = magnetization;
		}	
	}

	//============================
	// Save data triplets to file 
	//============================

	cnpy::npy_save(argv[3], data_triplets, {saved_data_size, 3}, "a");

	//======================
	// Deallocate resources 
	//======================

#ifdef RENDERING
	if (munmap(frame_buffer, fb_size) == -1)
	{
		fprintf(stderr, "[ISING-MODEL] Unable to unmap frame buffer from address space\n");
		exit(EXIT_FAILURE);
	}

	close(fb0_fd);
#endif

	return EXIT_SUCCESS;
}