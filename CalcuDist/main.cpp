#include "sort.h"

int main(int argc, char**argv)
{
	String fname = GHMIDir + (String)"HPH-144GHMIs.csv";
	vector<string> names;
	double ChiDist[N][N] = { 0.0 };

	Sort S(fname, ChiDist, names);

	while (1)
	{
		cout << "Please type in the image number that you wanna search(1~200):\n";
		int m;
		cin >> m;//µÚm¸öÍ¼Ïñ

		if (m > 200 || m < 1)
		{
			cout << "Input Error!" << endl;
			continue;
		}


		S.SortAndDisplay(ChiDist, m, names);

		/*String Dir = IMGDir;
		Mat Origin = imread(Dir + names[m - 1]);
		resize(Origin, Origin, Size(200, 200));
		namedWindow("Origin");
		imshow("Origin", Origin);*/

		cout << "Image " << m << " is " << names[m - 1] << endl;
		waitKey(0);
	}


	system("pause");
	return 0;
}