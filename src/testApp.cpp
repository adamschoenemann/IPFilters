#include "testApp.h"
#include <math.h>

void greyscaleFilter(ubyte*, int);

////////////////////////////////////////////////////////////////////////////////
ubyte clamp(int value){
	ubyte ret = value;
	if(value > 255){
		ret = 255;
	}
	else if(value < 0){
		ret = 0;
	}
	return ret;
}

////////////////////////////////////////////////////////////////////////////////
bool inBorder(int x, int y, int w, int h, int radius){
	if(x < radius || x + radius + 1 > w || y < radius || y + radius + 1 > h)
		return false;
	return true;
}


////////////////////////////////////////////////////////////////////////////////
int getIndex(int x, int y, int w){
	return (y*w+x) * 3;
}

////////////////////////////////////////////////////////////////////////////////
void bubbleSort(int *array,int length)
{
    int i,j;
    for(i=0;i<length;i++)
    {
        for(j=0;j<i;j++)
        {
            if(array[i]>array[j])
            {
                int temp=array[i]; //swap
                array[i]=array[j];
                array[j]=temp;
            }

        }

    }

}

////////////////////////////////////////////////////////////////////////////////
void toBinaryImage(ubyte* input, ubyte* output, int w, int h, ubyte t){
	for (int y = 0; y < h; y++) {
		for(int x = 0; x < w; x++){
			int i = (y*w+x)*3;
			int val = 0;
			for(int j = 0; j < 3; j++){
				val += input[i+j];
			}
			val /= 3;
			val = (val > t) ? 1 : 0;
			int b = y*w+x;
			output[b] = val;

		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void fromBinaryImage(ubyte* bin, ubyte* output, int w, int h){
	for (int y = 0; y < h; y++) {
		for(int x = 0; x < w; x++){
			int i = (y*w+x)*3;
			int b = y*w+x;
			for(int j = 0; j < 3; j++){

				output[i+j] = (bin[b] > 0) ? 255 : 0;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void prepGreyscale(int w, int h, ubyte* output){
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			int i = (y*w+x) * 3;
			greyscaleFilter(output, i);
		}
	}

}

//////////////////////////////// INVERT ////////////////////////////////////////
void invertFilter(ubyte* output, int i){
	output[i+0] = 255 - output[i];
	output[i+1] = 255 - output[i+1];
	output[i+2] = 255 - output[i+2];
}

////////////////////////////// BRIGHTNESS //////////////////////////////////////
void brightnessFilter(ubyte* output, int i, int brightness){
	output[i+0] = clamp(output[i] + brightness);
	output[i+1] = clamp(output[i+1] + brightness);
	output[i+2] = clamp(output[i+2] + brightness);
}

///////////////////////////// CONTRAST /////////////////////////////////////////
void contrastFilter(ubyte* output, int i, double contrast){
	int threshold = 127;
	for(int j = 0; j < 3; j++){
		output[i+j] = clamp(output[i+j]+(output[i+j]-threshold)*contrast);
	}
}

///////////////////////// AUTOCONTRAST /////////////////////////////////////////
struct AutoContrastFilter {

	double a;
	ubyte b;

	void setup(int w, int h, ubyte* output){
		ubyte min = output[0];
		ubyte max = output[0];

		for(int y = 0; y < h; y++){
			for (int x = 0; x < w; x++){
				int i = (y * w + x) * 3; // Index
				ubyte avg = (output[i]+output[i+1]+output[i+2])/3;
				if (avg > max){
					max = avg;
				}
				else if(avg < min){
					min = avg;
				}

			}

		}
		ubyte dif = max-min;

		a = 255.0/dif;
		b = min * -a;
	}

	void apply(ubyte* output, int i){
		for(int j = 0; j < 3; j++){
			output[i+j] = clamp(output[i+j]*a+b);
		}

	}
};

/////////////////////////// GREYSCALE //////////////////////////////////////////
void greyscaleFilter(ubyte* output, int i){
	ubyte avg = (output[i] + output[i+1] + output[i+2]) / 3;
	output[i] = avg;
	output[i+1] = avg;
	output[i+2] = avg;
}

//////////////////////////// GAMMA /////////////////////////////////////////////
void gammaFilter(ubyte* output, int i, double gamma){
	for (int j= 0; j < 3; j++){
		output[i+j] = clamp(255.0 * powf(output[i+j]/255.0, gamma));
	}
}

////////////////////////// ALPHA ///////////////////////////////////////////////
void alphaBlend(ubyte* output, ubyte*  blendPixels, int i, double a){
	for (int j=0; j < 3; j++){
		output[i+j] = a * output[i+j] + (1-a) * blendPixels[i+j];
	}
}

/////////////////////////////// THRESHOLD //////////////////////////////////////
void thresholdFilter(ubyte* output, int i, ubyte t){
	greyscaleFilter(output, i);
	ubyte temp;

	if(output[i] > t){
		temp = 255;
	}
	else{
		temp = 0;
	}
	for(int j = 0; j < 3; j++){
		output[i+j] = temp;
	}
}

//////////////////////////// NOISE /////////////////////////////////////////////
void saltPepperFilter(ubyte* output, int i, double thres){

	double randNum = ofRandom(0, 1);

	if(randNum > thres){

		ubyte val;
		double randNum2 = ofRandom(0, 1);
		if (randNum2 > 0.5) {
			val = 255;
		}
		else{
			val = 0;
		}
		for(int j = 0; j < 3; j++){
			output[i+j] = val;
		}

	}

}

////////////////////////////// HUE SHIFT ///////////////////////////////////////
void hueShift(ubyte* output, int i, int hue){
	ofColor col(output[i], output[i+1], output[i+2]);
	col.setHue(int(col.getHue() + hue) % 255);
	output[i] = col.r;
	output[i+1] = col.g;
	output[i+2] = col.b;
}

////////////////////////////// SATURATION //////////////////////////////////////
void saturationFilter(ubyte* output, int i, double sat){
	ofColor col(output[i], output[i+1], output[i+2]);
	col.setSaturation(col.getSaturation() * sat);
	output[i] = col.r;
	output[i+1] = col.g;
	output[i+2] = col.b;
}

///////////////////////////// SEPIA ////////////////////////////////////////////
void sepiaFilter(ubyte* output, int i){
	ubyte red;
	ubyte green;
	ubyte blue;

	red = output[i];
	green = output[i+1];
	blue = output[i+2];

	output[i] = clamp((red * .393) + (green *.769) + (blue * .189));
	output[i+1] = clamp((red * .349) + (green *.686) + (blue * .168));
	output[i+2] = clamp((red * .272) + (green *.534) + (blue * .131));
}

///////////////////////////// LOG MAP //////////////////////////////////////////
struct LogImageFilter {
	double c[3];
	void setup(ubyte* output, int w, int h){
		ubyte maxVal[3] = {0,0,0};
		for(int y = 0; y < h; y++){
			for(int x = 0; x < w; x++){
				int i;
				i = (w * y + x) * 3;
				for(int j = 0; j < 3; j++){
					if(maxVal[j] < output[i+j]){
						maxVal[j] = output[i+j];
					}
				}
			}
		}
		for(int j = 0; j< 3; j++){
			c[j] = 255/log(1+maxVal[j]);

		}
	}

	void apply(ubyte* output, int i){
		for(int j = 0; j < 3; j++){
			output[i+j] = c[j] * log(1+output[i+j]);

		}
	}

};

/////////////////////////// EXMP MAP ///////////////////////////////////////////
struct ExpImageFilter{

	double c[3];
	void setup(ubyte* output, int w, int h, double k){
		ubyte maxVal[3] = {0,0,0};
		for(int y = 0; y < h; y++){
			for(int x = 0; x < w; x++){
				int i;
				i = (w * y + x) * 3;
				for(int j = 0; j < 3; j++){
					if(maxVal[j] < output[i+j]){
						maxVal[j] = output[i+j];
					}
				}
			}
		}
		for(int j = 0; j < 3; j++){
			c[j] = 255.0 / (pow(k, double(maxVal[j])) - 1.0);

		}
	}

	void apply(ubyte* output, int i, double k){
		for(int j = 0; j < 3; j++){
			output[i+j] = c[j] * (pow(k, double(output[i+j])) - 1.0);
		}
	}

};

////////////////////////////// ACOS MAP ////////////////////////////////////////
void acosFilter(ubyte* output, int i){
	for(int j = 0; j < 3; j++){
		output[i+j] = acos(output[i+j] / 255.0) * output[i+j];

	}
}

////////////////////////////// COLOR DEPTH /////////////////////////////////////
void colDepth(ubyte* output, int i, int depth){
	for(int j = 0; j < 3; j++){
		int newMax = powf(2, depth)-1;
		double ratio = output[i+j]/255.0;
		int temp = ratio * newMax;
		output[i+j] = temp * 255/newMax;

	}
}

////////////////////////////// NORMALIZED RGB //////////////////////////////////
void normalizedRgb(ubyte* output, int i){

	int total = (output[i]+output[i+1]+output[i+2]);
	double r;
	double g;
	double I;

	if(total == 0){
		r = g = I = 0;
	}
	else{
		r = output[i]/(double)total;
		g = output[i+1]/(double)total;
		I = total/3.0;
	}

	ubyte R = r*total;
	ubyte G = g*total;
	ubyte B = (1.0-r-g)*total;

	output[i] = R;
	output[i+1] = G;
	output[i+2] = B;
}

///////////////////////////// PINK DETECTOR ////////////////////////////////////
void pinkDetector(ubyte* output,int i){
	int total = output[i] + output[i+1] + output[i+2];
	int temp = 0;
	if((output[i] > output[i+1] + 20 && output[i] > output[i+2] + 20) && total > 50){
		temp = 255;
	}
	for(int j = 0; j < 3; j++){
		output[i+j] = temp;
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////// NEIGHBORHOOD OPERATIONS /////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////// MEAN BLUR //////////////////////////////////
void meanBlur(ubyte* input, ubyte* output, int x, int y, int w, int h, int radius){

	if(inBorder(x, y, w, h, radius) == false) return;

	int size = radius * 2 + 1;
	int n = size * size;
	int means[3] = {0,0,0};

	for(int ky = -radius; ky < radius + 1; ky++){
		for(int kx = -radius; kx < radius + 1; kx++){

			int i = getIndex(x + kx, y + ky, w);
			for(int j = 0; j < 3; j++){
				means[j] += input[i+j];
			}
		}
	}
	int i = (y * w + x) * 3;
	for(int j = 0; j < 3; j++){
		output[i+j] = means[j] / n;
	}
}

////////////////////////////// MEAN DISC BLUR 7x7 //////////////////////////////
void meanDiscBlur(ubyte* input, ubyte* output, int x, int y, int w, int h, int radius){


	if(inBorder(x, y, w, h, radius) == false) return;

	int size = radius * 2 + 1;
	int* kernel = new int[size*size];

	for (int y = 0; y < size; y++) {
		for (int x = 0; x < size; x++) {
			int dx = x - radius;
			int dy = y - radius;
			int len = sqrt(dx*dx + dy*dy);
			int val = 0;
			if (len < radius) {
				val = 1;
			}
			kernel[y*size+x] = val;
		}
	}
	int n = 0;
	for(int c = 0; c < size * size; c++){
		n += kernel[c];
	}
	int means[3] = {0,0,0};

	for(int ky = -radius; ky < radius + 1; ky++){
		for(int kx = -radius; kx < radius + 1; kx++){
			//int i = ((y + ky) * w + (x + kx)) * 3;
			int i = getIndex(x + kx, y + ky, w);
			for(int j = 0; j < 3; j++){
				means[j] += input[i+j] * kernel[(ky+radius)*size+(kx+radius)];
			}
		}
	}
	int i = (y * w + x) * 3;
	for(int j = 0; j < 3; j++){
		output[i+j] = means[j] / n;
	}

	delete [] kernel;
}

///////////////////////////// MEDIAN ///////////////////////////////////////////
void median(ubyte* input, ubyte* output, int x, int y, int w, int h, int radius){

	if(inBorder(x, y, w, h, radius) == false) return;
	int size = radius * 2 + 1;
	int length = size * size;

	int* reds = new int[length];
	int* greens = new int[length];
	int* blues = new int[length];

	int k = 0;
	for (int ky = -radius ; ky <= radius ; ky++){
		for (int kx = -radius ; kx <= radius ; kx++){

			int i = getIndex(x + kx, y + ky, w);
			reds[k] = input[i];
			greens[k] = input[i+1];
			blues[k] = input[i+2];

			k++;
		}
	}
	bubbleSort(reds, length);
	bubbleSort(greens, length);
	bubbleSort(blues, length);
	int i = getIndex(x, y, w);
	output[i] = reds[length / 2];
	output[i+1] = greens[length / 2];
	output[i+2] = blues[length / 2];

	delete [] reds;
	delete [] greens;
	delete [] blues;
}

//////////////////////////////// GAUSSIAN BLUR 7x7 /////////////////////////////
void gaussianBlur(ubyte* input, ubyte* output, int x, int y, int w, int h){

	const int kw = 7;
	const int kh = 7;
	float kernel[kw][kh] =
	{
		0.00000067, 0.00002292, 0.00019117, 0.00038771, 0.00019117, 0.00002292, 0.00000067,
		0.00002292, 0.00078633, 0.00655965, 0.01330373, 0.00655965, 0.00078633, 0.00002292,
		0.00019117, 0.00655965, 0.05472157, 0.11098164, 0.05472157, 0.00655965, 0.00019117,
		0.00038771, 0.01330373, 0.11098164, 0.22508352, 0.11098164, 0.01330373, 0.00038771,
		0.00019117, 0.00655965, 0.05472157, 0.11098164, 0.05472157, 0.00655965, 0.00019117,
		0.00002292, 0.00078633, 0.00655965, 0.01330373, 0.00655965, 0.00078633, 0.00002292,
		0.00000067, 0.00002292, 0.00019117, 0.00038771, 0.00019117, 0.00002292, 0.00000067
	};

	int radius = 3;

	if(inBorder(x, y, w, h, radius) == false) return;

	int size = radius * 2 + 1;
	float sum = 0;
	for(int y = 0; y < kh; y++){
		for(int x = 0; x < kw; x++){
			sum += kernel[y][x];
		}
	}
	int totals[3] = {0,0,0};

	for(int ky = -radius; ky < radius + 1; ky++){
		for(int kx = -radius; kx < radius + 1; kx++){

			int i = getIndex(x + kx, y + ky, w);
			for(int j = 0; j < 3; j++){
				totals[j] += input[i+j] * kernel[ky+radius][kx+radius];
			}
		}
	}
	int i = (y * w + x) * 3;
	for(int j = 0; j < 3; j++){
		output[i+j] = totals[j] / sum;
	}
}

//////////////////////////// SOBEL EDGE DETECTION //////////////////////////////
void sobelFilter(ubyte* input, ubyte* output, int x, int y, int w, int h){
	int radius = 1;

	if(inBorder(x, y, w, h, radius) == false) return;
	int size = radius * 2 + 1;


	int n = 0;

	int means[3] = {0,0,0};

	for(int ky = -radius; ky < radius + 1; ky++){
		for(int kx = -radius; kx < radius + 1; kx++){
			int c = kx * ((radius * 2) - abs(ky));
			n++;
			int i = getIndex(x + kx, y + ky, w);
			for(int j = 0; j < 3; j++){
				means[j] += input[i+j] * c;
			}
		}
	}

	int i = (y * w + x) * 3;
	for(int j = 0; j < 3; j++){
		output[i+j] = clamp(means[j]);
	}

}

///////////////////////////// PREWITT EDGE DETECTION ///////////////////////////
void prewittFilter(ubyte* input, ubyte* output, int x, int y, int w, int h){
	int radius = 1;

	if(inBorder(x, y, w, h, radius) == false) return;
	int size = radius * 2 + 1;

	int n = 0;

	int means[3] = {0,0,0};

	for(int ky = -radius; ky < radius + 1; ky++){
		for(int kx = -radius; kx < radius + 1; kx++){
			int c = kx;
			n++;
			int i = getIndex(x + kx, y + ky, w);
			for(int j = 0; j < 3; j++){
				means[j] += input[i+j] * c;
			}
		}
	}

	int i = (y * w + x) * 3;
	for(int j = 0; j < 3; j++){
		output[i+j] = clamp(means[j]);
	}

}

///////////////////////////// MIN FILTER ///////////////////////////////////////
void minFilter(ubyte* input, ubyte* output, int x, int y, int w, int h){

	int radius = 1;
	if(inBorder(x, y, w, h, radius) == false) return;

	int size = radius * 2 + 1;
	int n = size * size;
	const int length = 9;
	int reds[length];
	int greens[length];
	int blues[length];

	int k = 0;
	for (int ky = -radius ; ky <= radius ; ky++){
		for (int kx = -radius ; kx <= radius ; kx++){

			int i = getIndex(x + kx, y + ky, w);
			reds[k] = input[i];
			greens[k] = input[i+1];
			blues[k] = input[i+2];

			k++;
		}
	}
	bubbleSort(reds, length);
	bubbleSort(greens, length);
	bubbleSort(blues, length);
	int i = getIndex(x, y, w);
	output[i] = reds[0];
	output[i+1] = greens[0];
	output[i+2] = blues[0];
}

//////////////////////////// MAX FILTER ////////////////////////////////////////
void maxFilter(ubyte* input, ubyte* output, int x, int y, int w, int h){

	int radius = 1;
	if(inBorder(x, y, w, h, radius) == false) return;

	int size = radius * 2 + 1;
	int n = size * size;
	const int length = 9;
	int reds[length];
	int greens[length];
	int blues[length];

	int k = 0;
	for (int ky = -radius ; ky <= radius ; ky++){
		for (int kx = -radius ; kx <= radius ; kx++){

			int i = getIndex(x + kx, y + ky, w);
			reds[k] = input[i];
			greens[k] = input[i+1];
			blues[k] = input[i+2];

			k++;
		}
	}
	bubbleSort(reds, length);
	bubbleSort(greens, length);
	bubbleSort(blues, length);
	int i = getIndex(x, y, w);
	output[i] = reds[length - 1];
	output[i+1] = greens[length - 1];
	output[i+2] = blues[length - 1];
}

///////////////////////////// EMBOSS ///////////////////////////////////////////
void emboss(ubyte* input, ubyte* output, int x, int y, int w, int h){

	int radius = 1;
	if(inBorder(x, y, w, h, radius) == false) return;

	const int kw = 3;
	const int kh = 3;
	float kernel[kw][kh] =
	{
		{-1, 0, -1},
		{0,  4, 0},
		{-1, 0, -1}
	};

	if (x < radius || x + radius + 1 > w || y < radius || y + radius +1 > h){
		return;
	}

	int size = radius * 2 + 1;
	float sum = 0;
	for(int y = 0; y < kh; y++){
		for(int x = 0; x < kw; x++){
			sum += kernel[y][x];
		}
	}
	int totals[3] = {0,0,0};

	for(int ky = -radius; ky < radius + 1; ky++){
		for(int kx = -radius; kx < radius + 1; kx++){

			int i = getIndex(x + kx, y + ky, w);
			for(int j = 0; j < 3; j++){
				totals[j] += input[i+j] * kernel[ky+radius][kx+radius];
			}
		}
	}
	int i = (y * w + x) * 3;
	for(int j = 0; j < 3; j++){
		output[i+j] = totals[j] + 127;
	}
}

////////////////////////////////////////////////////////////////////////////////
/////////////////////      MORPHOLOGY       ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////// ERODE /////////////////////////////////////////
void erode(ubyte* input, ubyte* output, int x, int y, int w, int h, int radius){

	if(inBorder(x, y, w, h, radius) == false) return;
	int size = radius * 2 + 1;

	bool* se = new bool[size*size];
	for(int i = 0; i < size*size; i++){
		se[i] = 1;
	}

	int k = 0;
	int j = y*w+x;
	for(int ky = -radius; ky < radius + 1; ky++){
		for(int kx = -radius; kx < radius + 1; kx++){
			int i = (y+ky)*w + (x+kx);

			if((input[i] & se[k]) != 1){
				output[j] = 0;
				delete [] se;
				return;
			}
			k++;
		}
	}

	delete [] se;

}

//////////////////////////// DILATE ////////////////////////////////////////////
void dilate(ubyte* input, ubyte* output, int x, int y, int w, int h, int radius){

	if(inBorder(x, y, w, h, radius) == false) return;

	int size = radius * 2 + 1;

	bool* se = new bool[size*size];
	for(int i = 0; i < size*size; i++){
		se[i] = 1;
	}

	int k = 0;
	int j = (y*w+x);
	for(int ky = -radius; ky < radius + 1; ky++){
		for(int kx = -radius; kx < radius + 1; kx++){
			int i = (y+ky)*w + (x+kx);

			if((input[i] & se[k]) == 1){
				output[j] = 1;
				delete [] se;
				return;
			}
			k++;
		}
	}

	delete [] se;
}

/////////////////////////////// CLOSE //////////////////////////////////////////
void close(ubyte* input, ubyte* output, int w, int h, int radius){
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			dilate(input, output, x, y, w, h, radius);
		}
	}
	memcpy(input, output, w*h);
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			erode(input, output, x, y, w, h, radius);
		}
	}
}

/////////////////////////// OPEN  //////////////////////////////////////////////
void open(ubyte* input, ubyte* output, int w, int h, int radius){
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			erode(input, output, x, y, w, h, radius);
		}
	}
	memcpy(input, output, w*h);
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			dilate(input, output, x, y, w, h, radius);
		}
	}
}

//////////////////////////// BOUNDARY //////////////////////////////////////////
void boundary(ubyte* input, ubyte* output, int w, int h, int radius){
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			erode(input, output, x, y, w, h, radius);
		}
	}

	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			int i = (y*w+x);
			output[i] = input[i] - output[i];
		}
	}
}
////////////////////////////////////////////////////////////////////////////////
//////////////////////////    BLOB ANALYSIS   //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// GRASS FIRE AND MORE IN HEADER FILE

/////////////////////////// BLUE COLOR TRACKING ///////////////////////////////
struct IsolateBlue {

	void trackColor(ubyte* input, ubyte* bin, int w, int h){
		for(int y = 0; y < h; y++){
			for(int x = 0; x < w; x++){
				int j = (y*w+x);
				int i = j * 3;
				ubyte* r = &input[i+0];
				ubyte* g = &input[i+1];
				ubyte* b = &input[i+2];
				double brRatio = double(*r) / *b;
				double grRatio = double(*r) / *g;
				double bgRatio = double(*g) / *b;
				unsigned int total = *r + *g + *b;
				if(brRatio < 0.3 && total > 40 && total < 600){
					bin[j] = 1;
				} else {
					bin[j] = 0;
				}
			}
		}
	}

};


/////////////////////////// BLOB ANALYSIS //////////////////////////////////////
BLOB analyzeBlob(ubyte* blob, int w, int h, int label){

	BLOB ret;
	ret.label = label;
	int prev = 0;
	for(int y = 0; y < h; y++){

		for(int x = 0; x < w; x++){
			int i = w*y+x;
			if(blob[i] == label){
				if(prev != label){
					ret.perim++;
				}
				prev = label;
				ret.centroid.x += x;
				ret.centroid.y += y;
				ret.box.set(x, y);
				ret.area++;
			} else {
				if(prev == label){
					ret.perim++;
					prev = 0;
				}
			}
		}

	}
	if(ret.area > 0){
		ret.centroid.x /= ret.area;
		ret.centroid.y /= ret.area;
		return ret;
	}

}

////////////////////////////////////////////////////////////////////////////////
/////////////////////    GEOMETRIC TRANSFORMATIONS  ////////////////////////////
////////////////////////////////////////////////////////////////////////////////


//////////////////////////// ROTATION //////////////////////////////////////////
void rotate(ubyte* input, ubyte* output, int w, int h, int rot){
	double rad = M_PI_2 * rot;
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){

			myPoint outPos(x, y);
			myPoint inPos(x, y);
			inPos.rotate(-rad, myPoint(w / 2, h / 2));
			if(inPos.x < 0 || inPos.x > w || inPos.y < 0 || inPos.y > h){
				continue;
			}
			int inIndex = (inPos.y * w + inPos.x) * 3;

			int outIndex = (y*w+x) * 3;
			for(int d = 0; d < 3; d++){
				output[outIndex + d] = input[inIndex + d];
			}
		}
	}
}

//////////////////////////// SCALE /////////////////////////////////////////////
void scale(ubyte* input, ubyte* output, int w, int h, double s){

	for(int y = 0; y < int(h * s); y++){
		for(int x = 0; x < int(w * s); x++){
			int i = (y*w+x) * 3;
			if(i > w*h*3) continue;
			int _x = x * 1.0/s;
			int _y = y * 1.0/s;
			int j = (_y*w+_x) * 3;
			for(int d = 0; d < 3; d++){
				output[i+d] = input[j+d];
			}
		}
	}
}

//////////////////////////// CROP //////////////////////////////////////////////
void crop(ubyte* input, ubyte* output, int w, int h, myPoint A, myPoint B){
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			if(ofInRange(x, A.x, B.x) && ofInRange(y, A.y, B.y)){
				int i = (y*w+x) * 3;
				for(int d = 0; d < 3; d++){
					output[i+d] = input[i+d];
				}
			}

		}
	}
}

////////////////////////////// TILE ////////////////////////////////////////////
void tile(ubyte* input, ubyte* output, int w, int h, int offsetX, int offsetY){
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){

			int in = (y*w+x) * 3;
			int nx = (x + offsetX);
			int ny = (y + offsetY);

			nx = (nx > w) ? nx - w : nx;
			ny = (ny > h) ? ny - h : ny;
			nx = (nx < 0) ? nx + w : nx;
			ny = (ny < 0) ? ny + h : ny;
			int out = (ny*w + nx) * 3;
			for(int d = 0; d < 3; d++){
				output[out+d] = input[in+d];
			}


		}
	}
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////     MAIN FUNCTIONS    //////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/////////////////////////// POINT PROCESS //////////////////////////////////////
void testApp::pointProcess(){

	pointInput = vidGrabber.getPixels();
	memcpy(pointOutput, pointInput, w*h*3);
//	AutoContrastFilter acf;
//	acf.setup(w, h, pointInput);
//	LogImageFilter lif;
//	lif.setup(pointOutput, w, h);
//	ExpImageFilter eif;
//	eif.setup(pointInput, w, h, 1.01);

	for(int y = 0; y < h; y++){
		for (int x = 0; x < w; x++){
			int i = (y * w + x) * 3; // Index

//				invertFilter(pointOutput, i);
//				brightnessFilter(pointOutput, i, 69);
//				contrastFilter(pointOutput, i, 0.6);
//				acf.apply(pointOutput, i);
//				greyscaleFilter(pointOutput, i);
//				gammaFilter(pointOutput, i, 0.45);
//				alphaBlend(pointOutput, blendPixels, i, 0.5);
//				thresholdFilter(pointOutput, i, 80);
//				saltPepperFilter(pointOutput, i, 0.9);
//				hueShift(pointOutput, i, 80);
//				saturationFilter(pointOutput, i, 2.5);
//				sepiaFilter(pointOutput, i);
//				lif.apply(pointOutput, i);
//				eif.apply(pointOutput, i, 1.01);
//				acosFilter(pointOutput, i);
				colDepth(pointOutput, i, 4);
//				normalizedRgb(pointOutput, i);
//				pinkDetector(pointOutput, i);


		}
	}
	videoTexture1.loadData(pointOutput, w, h, GL_RGB);
}

///////////////////////// NEIGHBORHOOD PROCESS /////////////////////////////////
void testApp::neighborhoodProcess(){

	pointInput = vidGrabber.getPixels();
	memcpy(pointOutput, pointInput, w*h*3);

	for(int y = 0; y < h; y++){
		for (int x = 0; x < w; x++){
			int i = (y * w + x) * 3; // Index
//			meanBlur(pointInput, pointOutput, x, y, w, h, 7);
//			median(pointInput, pointOutput, x, y, w, h, 2);
//			gaussianBlur(pointInput, pointOutput, x, y, w, h);
//			meanDiscBlur(pointInput, pointOutput, x, y, w, h, 7);
//			sobelFilter(pointInput, pointOutput, x, y, w, h);
//			prewittFilter(pointInput, pointOutput, x, y, w, h);
//			minFilter(pointInput, pointOutput, x, y, w, h);
//			maxFilter(pointInput, pointOutput, x, y, w, h);
			emboss(pointInput, pointOutput, x, y, w, h);
		}
	}

	videoTexture1.loadData(pointOutput, w, h, GL_RGB);
}

//////////////////////////// MORPHOLOGY ////////////////////////////////////////
void testApp::morphology(){
	pointInput = vidGrabber.getPixels();

	// Convert to Binary
	toBinaryImage(pointInput, binInput, w, h, 100);
	fromBinaryImage(binInput, blackPixels, w, h);
	memcpy(binOutput, binInput, w*h);

//	------------------------Not Compound----------------------------------------
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
//			erode(binInput, binOutput, x, y, w, h, 3);
//			dilate(binInput, binOutput, x, y, w, h, 3);
		}
	}

//------------------------------Compound----------------------------------------

//	open(binInput, binOutput, w, h, 3);

//	close(binInput, binOutput, w, h, 3);

	boundary(binInput, binOutput, w, h, 3);

	fromBinaryImage(binOutput, binDisplay, w, h);
	videoTexture1.loadData(blackPixels, w, h, GL_RGB);
	videoTexture2.loadData(binDisplay, w, h, GL_RGB);
}

/////////////////////////// BLOB ANALYSIS //////////////////////////////////////
void testApp::blobAnalysis(){

	pointInput = vidGrabber.getPixels();
	IsolateBlue iso;
	iso.trackColor(pointInput, binOutput, w, h);

	memcpy(binInput, binOutput, w*h);
//	memcpy(binOutput, binInput, w*h);

	close(binInput, binOutput, w, h, 4);
	memcpy(binInput, binOutput, w*h);

	// BLOB analysis -----------------------------------------------------------
	gf->startFire(binOutput, w, h);
	blob = analyzeBlob(gf->output, w, h, 1);
	cout << "centroid: " << blob.centroid.x << ", " << blob.centroid.y << endl;
	cout << "BB Ratio: " << blob.box.getRatio() << endl;
	cout << "Blob area: " << blob.area << endl;
	cout << "Blob perimeter: " << blob.perim << endl;
	cout << "Blob circularity: " << blob.getCircularity() << endl;
	cout << "Blob compactness: " << blob.getCompactness() << endl;

	GrassFire::blobToImage(gf->output, pointInput, w, h);
	videoTexture2.loadData(pointInput, w, h, GL_RGB);

}

////////////////////////// GEOMETRIC TRANSFORMATIONS ///////////////////////////
void testApp::transformations(){

	geomInput = vidGrabber.getPixels();
	memset(geomOutput, 0, w*h*3);

//	rotate(geomInput, geomOutput, w, h, 1);
//	scale(geomInput, geomOutput, w, h, 1.5);
//	crop(geomInput, geomOutput, w, h, myPoint(10, 10), myPoint(200, 100));
	tile(geomInput, geomOutput, w, h, -100, -100);

	geomTexture.loadData(geomOutput, w, h, GL_RGB);
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////// OPENFRAMEWORKS APP ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void testApp::setup(){

	w = 320;
	h = 240;

	image.loadImage("data/lars_reng.jpg");
	blendPixels = image.getPixels();

	vidGrabber.setVerbose(true); //Bool
	vidGrabber.initGrabber(w, h, true);

	videoTexture1.allocate(w,h, GL_RGB);
	videoTexture2.allocate(w,h, GL_RGB);
	geomTexture.allocate(w,h, GL_RGB);

	blackPixels = new ubyte[w*h*3];
	memset(blackPixels, 0, w*h*3);

	pointOutput = new ubyte[w*h*3];

	binInput = new ubyte[w*h];
	binOutput = new ubyte[w*h];
	binDisplay = new ubyte[w*h*3];

	geomInput = new ubyte[w*h*3];
	geomOutput = new ubyte[w*h*3];

	gf = new GrassFire(w, h);

	videoTexture1.loadData(blackPixels, w, h, GL_RGB);
	videoTexture2.loadData(blackPixels, w, h, GL_RGB);
	geomTexture.loadData(blackPixels, w, h, GL_RGB);

	ofBackground(255);
}

void testApp::update(){


	vidGrabber.grabFrame();

	if (vidGrabber.isFrameNew()){

//		pointProcess();
//		neighborhoodProcess();
//		morphology();
//		blobAnalysis();
		transformations();

	}
}

//--------------------------------------------------------------
void testApp::draw(){

	vidGrabber.draw(0, 0);
	videoTexture1.draw(w, 0);
	videoTexture2.draw(0, h);
	geomTexture.draw(w, h);

	if(gf->label > 0){
		ofPushStyle();
		ofCircle(blob.centroid.x, blob.centroid.y + h, 5);
		ofSetColor(255, 0, 0);
		ofNoFill();
		ofRect(blob.box.min.x, blob.box.min.y + h, blob.box.width(), blob.box.height());
		myPoint cBox = blob.box.getCenter();
		ofCircle(cBox.x, cBox.y + h, 5);
		ofPopStyle();
	}
	ofPushStyle();
	ofSetColor(0);
	ofDrawBitmapString(ofToString(ofGetFrameRate()), 30, 580);
	ofPopStyle();
}

//--------------------------------------------------------------
void testApp::keyPressed  (int key){

}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){

}
