// ========================================================================
// The Ising Model                                                         
// No Copyright. Vladislav Aleinik 2019                                    
// ========================================================================

#include <cstdlib>
#include <cstdio>
#include <cstring>

// ========================================================================
// Configuration File Work                                                
// ========================================================================

#include "Vector.hpp"

double temperature;
Vector externalField;
double interactivity;
double magnetic_moment;
size_t saved_data_samples;
size_t burn_in_samples;
size_t mc_iters_per_sample;
size_t iters_per_render_frame;

void read_config(const char* filename)
{
	FILE* conf_file = std::fopen(filename, "r");
	if (conf_file == NULL)
	{
		fprintf(stderr, "[ISING_MODEL] Unable to open config file\n");
		exit(EXIT_FAILURE);
	}

	if (fscanf(conf_file,            "temperature %lf\n",            &temperature) != 1)            temperature = 100.0;
	if (fscanf(conf_file,                  "field %lf\n",        &externalField.z) != 1)        externalField.z = 0.0;
	if (fscanf(conf_file,          "interactivity %lf\n",          &interactivity) != 1)          interactivity = 1.0;
	if (fscanf(conf_file,        "magnetic_moment %lf\n",        &magnetic_moment) != 1)        magnetic_moment = 1.0;
	if (fscanf(conf_file,     "saved_data_samples %zu\n",     &saved_data_samples) != 1)     saved_data_samples = 100;
	if (fscanf(conf_file,        "burn_in_samples %zu\n",        &burn_in_samples) != 1)        burn_in_samples = 20;
	if (fscanf(conf_file,    "mc_iters_per_sample %zu\n",    &mc_iters_per_sample) != 1)    mc_iters_per_sample = 100000;
	if (fscanf(conf_file, "iters_per_render_frame %zu\n", &iters_per_render_frame) != 1) iters_per_render_frame = 50000;

	interactivity *= 1.6e-19; // Joules

	externalField.x = 0.0;
	externalField.y = 0.0;

	fclose(conf_file);
}

#include "Model.hpp"

// ========================================================================
// Initial Spin States                                                     
// ========================================================================

long double getInteractivity(int x, int y)
{
	static std::random_device rd;
	static std::mt19937 gen{rd()};
	static std::uniform_real_distribution<double> distribution(-1.0, 1.0);

	return interactivity;
}

int getStateX(int x, int y)
{
	return rand() % STATE_GRAPH_SIZE_X;
}

int getStateY(int x, int y)
{
	return rand() % STATE_GRAPH_SIZE_Y;
}

// ========================================================================
#ifdef RENDERING
// ========================================================================

#include <fcntl.h>
#include <unistd.h>

#include <sys/mman.h>
#include <sys/ioctl.h>
#include <linux/fb.h>
#include <linux/kd.h>

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		fprintf(stderr, "[ISING-MODEL] Expected input: model <config file> <output file>\n");
		exit(EXIT_FAILURE);
	}
	
	//=========================
	// Read Configuration File 
	//=========================
	
	read_config(argv[1]);

	// Fixing units:

	externalField.z *= magnetic_moment;
	temperature     *= 1.38e-23;

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

	//==============
	// Calculations 
	//==============

	Lattice isingModel = Lattice(sizeX/8, sizeY/8, getInteractivity, getStateX, getStateY);

	while (true)
	{
		for (size_t i = 0; i < iters_per_render_frame; ++i)
		{
			isingModel.metropolisStep();
		}

		for (uint_fast16_t x = 0; x < sizeX/8; ++x)
		{
			for (uint_fast16_t y = 0; y < sizeY/8; ++y)
			{
				LatticePoint curPoint = isingModel.get(x, y);
				Vector spin = isingModel.stateGraph.get(curPoint.stateX, curPoint.stateY);

				for (size_t dx = 0; dx < 8; ++dx) {
				for (size_t dy = 0; dy < 8; ++dy) {
					uint_fast16_t pix_offset = (8*x+dx) * bytes_per_pixel +
					                           (8*y+dy) * bytes_per_line;
					frame_buffer[pix_offset + offset_red  ] = (1.0 + spin.x) * 127;;
					frame_buffer[pix_offset + offset_green] = (1.0 + spin.z) * 127;
					frame_buffer[pix_offset + offset_blue ] = (1.0 + spin.y) * 127;;
				}}
			}
		}
	}

	//======================
	// Deallocate resources 
	//======================

	if (munmap(frame_buffer, fb_size) == -1)
	{
		fprintf(stderr, "[ISING-MODEL] Expected input: model <config file> <output file>\n");
		exit(EXIT_FAILURE);
	}

	close(fb0_fd);

	return EXIT_SUCCESS;
}

// ========================================================================
#else
// ========================================================================

#include <thread>
#include <functional>
#include <vector>

#include "vendor/cnpy/cnpy.h"

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		fprintf(stderr, "[ISING-MODEL] Expected input: model <temperature> <field> <filename>\n");
		exit(EXIT_FAILURE);
	}
	
	//=========================
	// Read Configuration File 
	//=========================

	read_config(argv[1]);
	
	// Fixing units:

	double oldFieldZ = externalField.z;
	double oldT      = temperature;

	externalField.z *= magnetic_moment;
	temperature *= 1.38e-23;

	//===================================
	// Allocate buffer for data triplets 
	//===================================

	double* data_points = (double*) calloc(2 * saved_data_samples, sizeof(*data_points));
	if (data_points == NULL)
	{
		fprintf(stderr, "[ISING-MODEL] Unable to allocate array for saved data!\n");
		exit(EXIT_FAILURE);
	}

	//==============
	// Calculations 
	//==============

	printf("Computing for T=%lf H=%lf\n", oldT, oldFieldZ);

	int sizeX = 200;
	int sizeY = 200;
	Lattice isingModel = Lattice(sizeX, sizeY, getInteractivity, getStateX, getStateY);

	for (size_t iteration = 0, cur_saved_data = 0; iteration < burn_in_samples + saved_data_samples; ++iteration)
	{
		printf("\rComputation in progress: %02.0f%%", 100.0 * cur_saved_data/saved_data_samples);
		fflush(stdout);

		for (size_t iter = 0; iter < mc_iters_per_sample; ++iter)
			isingModel.metropolisStep();

		if (burn_in_samples <= iteration && cur_saved_data < saved_data_samples)
		{
			data_points[2 * cur_saved_data + 0] = magnetic_moment * isingModel.calculateMagnetization();
			data_points[2 * cur_saved_data + 1] = isingModel.calculateEnergy();

			++cur_saved_data;
		}
	}
	
	printf("\rComputation in progress: %02.0f%%", 100.0);
	printf("\nComputation completed!\n");

	//============================
	// Save data triplets to file 
	//============================

	cnpy::npy_save(argv[2], data_points, {saved_data_samples, 2}, "w");
	
	delete[] data_points;

	return EXIT_SUCCESS;
}

#endif