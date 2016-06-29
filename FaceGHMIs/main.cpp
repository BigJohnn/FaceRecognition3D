//Copyright(c) 2016, Hou Peihong
//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met :
//
//1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
//2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation 
//and / or other materials provided with the distribution.
//
//3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software
//without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
//BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#include <stdio.h>
#include <stdlib.h>
#include <iostream>             // for std::cout 
#include <string>
#include <fstream>
#include <algorithm>
#include <pcl/io/pcd_io.h>      // header that contains the definitions for PCD I/O operations 
#include <pcl/point_types.h>    // header that contains definitions for several PointT type structures 
#include <boost/filesystem.hpp> // includes all needed Boost.Filesystem declarations to find all data base vrml files 
#include <iterator>
#include <vector>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/conditional_removal.h>
#include <math.h>
#include <boost/filesystem.hpp> 

using namespace std;
namespace fs = boost::filesystem;

#define N 3
#define imagenumMAX 3500 //电脑异，此处最值亦不同矣.
#define InvDim 10
#define _STR(s) #s
#define STR(s) _STR(s)
#define normpara 100 

typedef double Invtype;

int get_filenames(const std::string& dir, std::vector<std::string>& filenames)
{
	fs::path path(dir);
	if (!fs::exists(path))
	{
		return -1;
	}

	fs::directory_iterator end_iter;
	for (fs::directory_iterator iter(path); iter != end_iter; ++iter)
	{
		if (fs::is_regular_file(iter->status()))
		{
			filenames.push_back(iter->path().string());
		}

		if (fs::is_directory(iter->status()))
		{
			get_filenames(iter->path().string(), filenames);
		}
	}

	return filenames.size();
}

//Calculate Hermite Polynomials by its recursive relationship.
double H(int p, double x)
{
	if (p == 0)
		return 1;
	else if (p == 1)
		return 2 * x;
	else
		return 2 * x*H(p - 1, x) - 2 * (p - 1)*H(p - 2, x);
}

//Function to split a string 
//from source: http://shadow2531.com/cpp/splitter.cpp
/////////////////////////////////////////////////////// 
inline vector<string> split(const string& s, const string& f)
{
	vector<string> temp;

	//if no seperator than return string 
	if (f.empty()) {
		temp.push_back(s);
		return temp;
	}

	typedef string::const_iterator iter;
	const iter::difference_type f_size(distance(f.begin(), f.end()));
	iter i(s.begin());

	for (iter pos; (pos = search(i, s.end(), f.begin(), f.end())) != s.end();)
	{
		temp.push_back(string(i, pos));
		advance(pos, f_size);
		i = pos;
	}

	temp.push_back(string(i, s.end()));

	return temp;
}


vector<Invtype> Calculate18GHMIs4SingleImage(vector<Invtype> &Invariant, pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud, double sigma)
{
	double M[4][4][4] = { 0.0 };
	double xx, yy,zz;
	double a;
	//Calculate--GHMs--000~333
	for (int r = 0; r <= N; r++)
	{
		for (int q = 0; q <= N; q++)
		{
			for (int p = 0; p <= N; p++)
			{
				for (size_t i = 0; i < cloud->size(); i++)
				{
					xx = cloud->points[i].x*1.0 / normpara;
					yy = cloud->points[i].y*1.0 / normpara;
					zz = cloud->points[i].z*1.0 / normpara;

					a = cloud->points[i].rgb;
					M[p][q][r] += cloud->points[i].rgb * pow(10, 39) / 255 * H(p, xx / sigma)*H(q, yy / sigma)*H(r, zz / sigma)*exp(-(xx*xx + yy*yy + zz*zz) / 2 / (sigma*sigma));
				}
			}
		}
	}

	Invtype Inv[InvDim] = { 0.0 };

	double u000 = M[0][0][0];
	double u200 = M[2][0][0];
	double u020 = M[0][2][0];
	double u002 = M[0][0][2];
	double u110 = M[1][1][0];
	double u101 = M[1][0][1];
	double u011 = M[0][1][1];
	double u300 = M[3][0][0];
	double u030 = M[0][3][0];
	double u003 = M[0][0][3];
	double u210 = M[2][1][0];
	double u201 = M[2][0][1];
	double u120 = M[1][2][0];
	double u102 = M[1][0][2];
	double u021 = M[0][2][1];
	double u012 = M[0][1][2];
	double u111 = M[1][1][1];

	//I1~I5
	Inv[0] = (u200 + u020 + u002) / u000 / u000; 

	Inv[1] = (u200*u200 + u020*u020 + u002*u002 + 2 * u110*u110 + 2 * u101*u101 + 2 * u011*u011
		) / u000 / u000 / u000 / u000;

	Inv[2] = (u200*u200*u200 + 3 * u200*u110*u110 + 3 * u200*u101*u101
		+ 3 * u110*u110*u020 + 3 * u101*u101*u002 + u020*u020*u020 + 3 * u020*u011*u011
		+ 3 * u011*u011*u002 + u002*u002*u002 + 6 * u110*u101*u011) / u000 / u000 / u000 / u000 / u000 / u000;

	Inv[3] = (u300*u300 + u030*u030 + u003*u003 + 3 * u210*u210 + 3 * u201*u201 + 3 * u120*u120 + 3 * u102*u102 + 3 * u021*u021 + 3 * u012*u012
		+ 6 * u111*u111) / u000 / u000 / u000 / u000 / u000;

	Inv[4] = (u300*u300 + 2 * u300*u120 + 2 * u300*u102 + 2 * u210*u030 + 2 * u201*u003 + u030*u030 + 2 * u030*u012
		+ 2 * u021*u003 + u003*u003 + u210*u210 + 2 * u210*u012 + u201*u201 + 2 * u201*u021 + u120*u120 + 2 * u120*u102
		+ u102*u102 + u021*u021 + u012*u012) / u000 / u000 / u000 / u000 / u000;

	//I9
	Inv[5] = (u200*u300*u300 + 2 * u110*u300*u210 + 2 * u110*u120*u030 + 2 * u101*u300*u201 + 2 * u101*u102*u003
		+ u020*u030*u030 + 2 * u011*u030*u021 + 2 * u011*u012*u003 + u002*u003*u003 + 2 * u200*u210*u210
		+ 2 * u200*u201*u201 + u200*u120*u120 + 2 * u200*u111*u111 + u200*u102*u102 + 4 * u110*u210*u120
		+ 4 * u110*u201*u111 + 4 * u110*u111*u021 + 2 * u110*u102*u012 + 4 * u101*u210*u111
		+ 4 * u101*u201*u102 + 2 * u101*u120*u021 + 4 * u101*u111*u012 + u020*u210*u210 + 2 * u020*u120*u120
		+ 2 * u020*u111*u111 + 2 * u020*u021*u021 + u020*u012*u012 + 2 * u011*u210*u201 + 4 * u011*u120*u111
		+ 4 * u011*u111*u102 + 4 * u011*u021*u012 + u002*u201*u201 + 2 * u002*u111*u111 + 2 * u002*u102*u102
		+ u002*u021*u021 + 2 * u002*u012*u012) / u000 / u000 / u000 / u000 / u000 / u000 / u000;

	//I10
	Inv[6] = (u200*u300*u300 + u200*u300*u120 + u200*u300*u102 + u200*u210*u030 + u200*u201*u003
		+ 2 * u110*u300*u210 + 2 * u110*u120*u030 + 2 * u110*u111*u003 + 2 * u101*u300*u201
		+ 2 * u101*u111*u030 + 2 * u101*u102*u003 + u020*u300*u120 + u020*u210*u030 + u020*u030*u030
		+ u020*u030*u012 + u020*u021*u003 + 2 * u011*u300*u111 + 2 * u011*u030*u021
		+ 2 * u011*u012*u003 + u002*u300*u102 + u002*u201*u003 + u002*u030*u012 + u002*u021*u003
		+ u002*u003*u003 + u200*u210*u210 + u200*u210*u012 + u200*u201*u201 + u200*u201*u021 + 4 * u110*u210*u120
		+ 2 * u110*u210*u102 + 2 * u110*u201*u111 + 2 * u110*u120*u012 + 2 * u110*u111*u021
		+ 2 * u101*u210*u111 + 2 * u101*u201*u120 + 4 * u101*u201*u102 + 2 * u101*u111*u012
		+ 2 * u101*u102*u021 + u020*u201*u021 + u020*u120*u120 + u020*u120*u102 + u020*u021*u021
		+ 2 * u011*u210*u021 + 2 * u011*u201*u012 + 2 * u011*u120*u111 + 2 * u011*u111*u102
		+ 4 * u011*u021*u012 + u002*u210*u012 + u002*u120*u102 + u002*u102*u102 + u002*u012*u012) / u000 / u000 / u000 / u000 / u000 / u000 / u000;

	//I11
	Inv[7] = (u200*u300*u300 + 2 * u200*u300*u120 + 2 * u200*u300*u102 + 2 * u110*u300*u210 + 2 * u110*u300*u030
		+ 2 * u110*u300*u012 + 2 * u110*u120*u030 + 2 * u110*u102*u030 + 2 * u101*u300*u201
		+ 2 * u101*u300*u021 + 2 * u101*u300*u003 + 2 * u101*u120*u003 + 2 * u101*u102*u003
		+ 2 * u020*u210*u030 + u020*u030*u030 + 2 * u020*u030*u012 + 2 * u011*u210*u003 + 2 * u011*u201*u030
		+ 2 * u011*u030*u021 + 2 * u011*u030*u003 + 2 * u011*u012*u003 + 2 * u002*u201*u003
		+ 2 * u002*u021*u003 + u002*u003*u003 + u200*u120*u120 + 2 * u200*u120*u102 + u200*u102*u102
		+ 2 * u110*u210*u120 + 2 * u110*u210*u102 + 2 * u110*u120*u012 + 2 * u110*u102*u012
		+ 2 * u101*u201*u120 + 2 * u101*u201*u102 + 2 * u101*u120*u021 + 2 * u101*u102*u021 + u020*u210*u210
		+ 2 * u020*u210*u012 + u020*u012*u012 + 2 * u011*u210*u201 + 2 * u011*u210*u021 + 2 * u011*u201*u012
		+ 2 * u011*u021*u012 + u002*u201*u201 + 2 * u002*u201*u021 + u002*u021*u021) / u000 / u000 / u000 / u000 / u000 / u000 / u000;

	//I12
	Inv[8] = (u200*u200*u300*u300 + 4 * u200*u110*u300*u210 + 4 * u200*u101*u300*u201 + 2 * u200*u020*u300*u120
		+ 2 * u200*u020*u210*u030 + 4 * u200*u011*u300*u111 + 2 * u200*u002*u300*u102
		+ 2 * u200*u002*u201*u003 + 4 * u110*u020*u120*u030 + 4 * u110*u002*u111*u003
		+ 4 * u101*u020*u111*u030 + 4 * u101*u002*u102*u003 + u020*u020*u030*u030 + 4 * u020*u011*u030*u021
		+ 2 * u020*u002*u030*u012 + 2 * u020*u002*u021*u003 + 4 * u011*u002*u012*u003 + u002*u002*u003*u003
		+ u200*u200*u210*u210 + u200*u200*u201*u201 + 4 * u200*u110*u210*u120 + 4 * u200*u110*u201*u111
		+ 4 * u200*u101*u210*u111 + 4 * u200*u101*u201*u102 + 2 * u200*u020*u201*u021
		+ 4 * u200*u011*u210*u021 + 4 * u200*u011*u201*u012 + 2 * u200*u002*u210*u012 + 4 * u110*u110*u210*u210
		+ 4 * u110*u110*u120*u120 + 8 * u110*u101*u210*u201 + 8 * u110*u101*u120*u111 + 8 * u110*u101*u111*u102
		+ 4 * u110*u020*u210*u120 + 4 * u110*u020*u111*u021 + 8 * u110*u011*u210*u111
		+ 8 * u110*u011*u120*u021 + 8 * u110*u011*u111*u012 + 4 * u110*u002*u210*u102
		+ 4 * u110*u002*u120*u012 + 4 * u101*u101*u201*u201 + 4 * u101*u101*u102*u102 + 4 * u101*u020*u201*u120
		+ 4 * u101*u020*u102*u021 + 8 * u101*u011*u201*u111 + 8 * u101*u011*u111*u021
		+ 8 * u101*u011*u102*u012 + 4 * u101*u002*u201*u102 + 4 * u101*u002*u111*u012 + u020*u020*u120*u120
		+ u020*u020*u021*u021 + 4 * u020*u011*u120*u111 + 4 * u020*u011*u021*u012 + 2 * u020*u002*u120*u102
		+ 4 * u011*u011*u021*u021 + 4 * u011*u011*u012*u012 + 4 * u011*u002*u111*u102 + 4 * u011*u002*u021*u012 + u002*u002*u102*u102
		+ u002*u002*u012*u012 + 4 * u110*u110*u111*u111 + 4 * u101*u101*u111*u111 + 4 * u011*u011*u111*u111) / u000 / u000 / u000 / u000 / u000 / u000 / u000 / u000 / u000;

	//I13
	Inv[9] = (u200*u200*u300*u300 + 2 * u200*u110*u300*u210 + 2 * u200*u110*u120*u030 + 2 * u200*u101*u300*u201
		+ 2 * u200*u101*u102*u003 + u110*u110*u300*u300 + u110*u110*u030*u030 + 2 * u110*u101*u030*u021
		+ 2 * u110*u101*u012*u003 + 2 * u110*u020*u300*u210 + 2 * u110*u020*u120*u030
		+ 2 * u110*u011*u300*u201 + 2 * u110*u011*u102*u003 + u101*u101*u300*u300 + u101*u101*u003*u003
		+ 2 * u101*u011*u300*u210 + 2 * u101*u011*u120*u030 + 2 * u101*u002*u300*u201
		+ 2 * u101*u002*u102*u003 + u020*u020*u030*u030 + 2 * u020*u011*u030*u021 + 2 * u020*u011*u012*u003
		+ u011*u011*u030*u030 + u011*u011*u003*u003 + 2 * u011*u002*u030*u021 + 2 * u011*u002*u012*u003 + u002*u002*u003*u003
		+ 2 * u200*u200*u210*u210 + 2 * u200*u200*u201*u201 + u200*u200*u120*u120 + 2 * u200*u200*u111*u111 + u200*u200*u102*u102
		+ 4 * u200*u110*u210*u120 + 4 * u200*u110*u201*u111 + 4 * u200*u110*u111*u021
		+ 2 * u200*u110*u102*u012 + 4 * u200*u101*u210*u111 + 4 * u200*u101*u201*u102
		+ 2 * u200*u101*u120*u021 + 4 * u200*u101*u111*u012 + 3 * u110*u110*u210*u210 + 2 * u110*u110*u201*u201
		+ 3 * u110*u110*u120*u120 + u110*u110*u102*u102 + 2 * u110*u110*u021*u021 + u110*u110*u012*u012 + 2 * u110*u101*u210*u201
		+ 4 * u110*u101*u120*u111 + 4 * u110*u101*u111*u102 + 4 * u110*u101*u021*u012
		+ 4 * u110*u020*u210*u120 + 4 * u110*u020*u201*u111 + 4 * u110*u020*u111*u021
		+ 2 * u110*u020*u102*u012 + 4 * u110*u011*u210*u111 + 4 * u110*u011*u201*u102
		+ 2 * u110*u011*u120*u021 + 4 * u110*u011*u111*u012 + 2 * u101*u101*u210*u210 + 3 * u101*u101*u201*u201 + u101*u101*u120*u120
		+ 3 * u101*u101*u102*u102 + u101*u101*u021*u021 + 2 * u101*u101*u012*u012 + 4 * u101*u011*u210*u120 + 4 * u101*u011*u201*u111
		+ 4 * u101*u011*u111*u021 + 2 * u101*u011*u102*u012 + 4 * u101*u002*u210*u111
		+ 4 * u101*u002*u201*u102 + 2 * u101*u002*u120*u021 + 4 * u101*u002*u111*u012 + u020*u020*u210*u210
		+ 2 * u020*u020*u120*u120 + 2 * u020*u020*u111*u111 + 2 * u020*u020*u021*u021 + u020*u020*u012*u012 + 2 * u020*u011*u210*u201
		+ 4 * u020*u011*u120*u111 + 4 * u020*u011*u111*u102 + 4 * u020*u011*u021*u012 + u011*u011*u210*u210
		+ u011*u011*u201*u201 + 2 * u011*u011*u120*u120 + 2 * u011*u011*u102*u102 + 3 * u011*u011*u021*u021 + 3 * u011*u011*u012*u012
		+ 2 * u011*u002*u210*u201 + 4 * u011*u002*u120*u111 + 4 * u011*u002*u111*u102
		+ 4 * u011*u002*u021*u012 + u002*u002*u201*u201 + 2 * u002*u002*u111*u111 + 2 * u002*u002*u102*u102 + u002*u002*u021*u021
		+ 2 * u002*u002*u012*u012 + 4 * u110*u110*u111*u111 + 4 * u101*u101*u111*u111 + 4 * u011*u011*u111*u111) / u000 / u000 / u000 / u000 / u000 / u000 / u000 / u000 / u000;

	double S = 0.0;

	for (int i = 0; i < InvDim; i++)
	{
		S += Inv[i] * Inv[i];
	}

	double D = sqrt(S);

	if (0 == D)
		D = 1;

	for (int i = 0; i < InvDim; i++)
	{
		Inv[i] /= D;
		cout << Inv[i] << endl;
		Invariant.push_back(Inv[i]);
	}
	return Invariant;
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr readWrl(string &file_name)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>); //initialize point cloud 
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZRGB>);

	fstream file;
	string s;
	string filename = file_name;

	//open data base file	   
	file.open(filename.c_str(), ios::in);

	//while not end of file 
	while (!file.eof())
	{
		//get line and save as string 
		getline(file, s);

		//if string = point [ (start of the coordinates in wrl file) 
		if (s.find("point [", 0) != string::npos)
		{
			//continue reading 
			getline(file, s);
			int j = 1; //int for cloud.width for unstructured point cloud 
			size_t i = 0; //point cloud element 

			//until string = ] (end of the coordinates) 
			while (s.find("]", 0) == string::npos)
			{
				//split string with seperator 
				vector<string> temp = split(s, " ");

				// Fill in the cloud data for unstructured point clouds 
				cloud->width = j;
				cloud->height = 1;
				cloud->is_dense = false;
				cloud->points.resize(cloud->width * cloud->height);

				for (; i < cloud->points.size(); ++i)
				{
					//connect cartesian coordinates from splitted vector 
					cloud->points[i].x = atof(temp[1].c_str());
					cloud->points[i].y = atof(temp[2].c_str());
					cloud->points[i].z = atof(temp[3].c_str()) + 1600;
				}

				//read next line 
				getline(file, s);
				j++;
			}
		}

		if (s.find("color [", 0) != string::npos)
		{
			//continue reading 
			getline(file, s);
			size_t i = 0; //point cloud element 

			//until string = ] (end of the coordinates) 
			for (; i < cloud->points.size(); ++i)
			{
				//split string with seperator 
				vector<string> temp = split(s, " ");

				if (s.find("]", 0) == string::npos)
				{
					cloud->points[i].r = (uint8_t)(atof(temp[1].c_str()) * 255);
					cloud->points[i].g = (uint8_t)(atof(temp[2].c_str()) * 255);
					cloud->points[i].b = (uint8_t)(atof(temp[3].c_str()) * 255);

					getline(file, s);
				}
				else break;
			}
			break;
		}
	}

	//close file 
	file.close();

	/*string sName = dir_iter->path().stem().string() + ".pcd";
	string sAll = sOut + "\\" + sName;
	pcl::io::savePCDFileASCII(sAll, cloud);
	std::cerr << endl << "Saved " << cloud.points.size() << " data points to " << sName;*/

	//Preprocessing
	pcl::RadiusOutlierRemoval<pcl::PointXYZRGB> outrem;
	// build the filter
	outrem.setInputCloud(cloud);
	outrem.setRadiusSearch(10);// 搜索半径
	outrem.setMinNeighborsInRadius(65); //有k个邻接点则保留
	// apply filter
	outrem.filter(*cloud_filtered);

	return cloud_filtered;
}

vector<Invtype> CalImgGHMIs(string filename)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud = readWrl(filename);

	vector<Invtype> Invariant;

	if (cloud->size() == 0)
	{
		Invariant.clear();
		return Invariant;
	}

	//Get Invariant whose dimension is 4*InvDim(18)*3≈216.
	double sigma = 0.2;
	for (int i = 0; i < 8; i++)
	{
		Invariant = Calculate18GHMIs4SingleImage(Invariant, cloud, sigma);
		sigma += 0.2;
	}

	/*this->feature= new double[];

	for (int i = 0; i<this->m_imgCount; i++)
	feature[i] = Invariant[i];*/

	return Invariant;
}


boost::shared_ptr<pcl::visualization::PCLVisualizer> simpleVis(pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud)
{
	// --------------------------------------------
	// -----Open 3D viewer and add point cloud-----
	// --------------------------------------------
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addPointCloud<pcl::PointXYZ>(cloud, "sample cloud");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "sample cloud");
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();
	return (viewer);
}

boost::shared_ptr<pcl::visualization::PCLVisualizer> rgbVis(pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud)
{
	// --------------------------------------------
	// -----Open 3D viewer and add point cloud-----
	// --------------------------------------------
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	viewer->setBackgroundColor(0, 0, 0);
	pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(cloud);
	viewer->addPointCloud<pcl::PointXYZRGB>(cloud, rgb, "sample cloud");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud");
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();
	return (viewer);
}



int main(int argc, char* argv[])
{
	char* suffix = ".csv";
	string Dir = (string)"HPH-" + "64GHMIs" + (string)suffix;
	ofstream fout(Dir.c_str());

	//searching in data base folder for all VRML object models to transform them later into pcl poinclouds 
	//using boost library 
	
	string sIn = argv[1];
	if (!argv[1]) return -1;

	cout << "Input directory: " << sIn << endl;
	
	vector<Invtype> Invariant;
	vector<string> filenames;

	int img_count = get_filenames(sIn.c_str(), filenames);
	
	vector< vector<Invtype> > Invs;
	vector<Invtype> inv;
	for (int i = 0; i < img_count; i++)  
	{
		inv = CalImgGHMIs(filenames[i]);

		if (inv.size() == 0) continue;

		Invs.push_back(inv);

		fout << filenames[i];

		for (size_t i = 0; i != inv.size(); i++)
			fout << "," << inv[i];
		fout << endl;
	}

	//boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;

	//viewer = rgbVis(cloud_filtered);


	////pcl::io::savePCDFileASCII("result.pcd", result);

	//while (!viewer->wasStopped())
	//{
	//	viewer->spinOnce(100);
	//	boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	//}

	system("pause");

	return 0;
}

