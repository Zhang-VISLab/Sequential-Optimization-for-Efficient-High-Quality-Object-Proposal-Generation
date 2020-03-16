#include "stdafx.h"
#include "DataSetVOC.h"
#include "FelzenSegment/segment-image.h"

inline int clamp(int v, int a, int b) { return (v<a ? a : v>b ? b : v); }
#define fast_max(x,y) (x - ((x - y) & ((x - y) >> (sizeof(int) * CHAR_BIT - 1))))
#define fast_min(x,y) (y + ((x - y) & ((x - y) >> (sizeof(int) * CHAR_BIT - 1))))

void boxesNms(ValStructVec<float, Vec4i> &boxes, float thr, int maxBoxes)
{
	boxes.sort();
	if (thr > .99)
		return;
	const int nBin = 60;
	const float step = 1 / thr, lstep = log(step);
	vector<ValStructVec<float, Vec4i> > kept(nBin + 1);
	int i = 0, j, k, n = (int)boxes.size(), m = 0, b;
	while (i < n && m < maxBoxes) {
		b = (boxes[i][2] - boxes[i][0] + 1) * (boxes[i][3] - boxes[i][1] + 1);
		bool keep = true;
		b = clamp(int(ceil(log(float(b)) / lstep)), 1, nBin);
		for (j = b - 1; j <= b + 1; j++) {
			int num = kept[j].size();
			for (k = 0; k < num; k++) {
				keep = DataSetVOC::interUnio(boxes[i], kept[j][k]) <= thr;
				if (!keep)
					goto LL;
			}
		}
LL:		if (keep) {
			kept[b].pushBack(boxes(i),boxes[i]);
			m++;
		}
		i++;
	}
	boxes.reserve(m);
	for (j = 1; j <= nBin; j++)
	    for (k = 0; k < kept[j].size(); k++)
			boxes.pushBack(kept[j](k), kept[j][k]);
    boxes.sort();
}

//SegmentRecursiveBox
void SegmentRecursiveBox(vector<Vec4i> &sp_boxes, ValStructVec<float, Vec4i> &init_boxes, const vecF thetas, const float beta, ValStructVec<float, Vec4i> &out_boxes)
{
	const int numS = sp_boxes.size();
	const int numI = init_boxes.size();
	const int numT = thetas.size();

    // preallocate memory    
    int *areaS = new int[numS];
	vector<pair<float, int> > overlap(numS);
    
    // some variables
    int i, j, k;
    int inner_count, inter_count;
    int x1, x2, y1, y2, mx1, mx2, my1, my2;
	int ox1, ox2, oy1, oy2, rx1, rx2, ry1, ry2;
    int areaI, area_inter, area_inner, areaM;
    float overlap_max, overlap_new, inter;
    int nnzbox = 0;        
    
    // box initialization
	for (i = 0; i < numI; i++)
		init_boxes(i) = numI - i;
    
    // areas of sp_boxes    
    for(i = 0; i < numS; i++)
		areaS[i] = (sp_boxes[i][2] - sp_boxes[i][0] + 1)*(sp_boxes[i][3] - sp_boxes[i][1] + 1);
    
    // stage 1: box alignment
	for (i = 0; i < numI; i++) {
		x1 = init_boxes[i][0], y1 = init_boxes[i][1];
		x2 = init_boxes[i][2], y2 = init_boxes[i][3];
        areaI = (x2 - x1 + 1) * (y2 - y1 + 1);
        rx1 = x1; ry1 = y1; rx2 = x2; ry2 = y2;
        
		inter_count = inner_count = 0;
        for (j = 0; j < numS; j++) {
            // overlap
			ox1 = fast_max(x1, sp_boxes[j][0]);
			ox2 = fast_min(x2, sp_boxes[j][2]);
			if (ox1 > ox2) continue;
			oy1 = fast_max(y1, sp_boxes[j][1]);
			oy2 = fast_min(y2, sp_boxes[j][3]);
			if (oy1 > oy2) continue;
            
            // inner regions
			if (ox1 == sp_boxes[j][0] && oy1 == sp_boxes[j][1] && ox2 == sp_boxes[j][2] && oy2 == sp_boxes[j][3]) {
				inner_count++;
                // combine regions
                if (inner_count == 1)
					rx1 = ox1, ry1 = oy1, rx2 = ox2, ry2 = oy2;
                else {
                    rx1 = fast_min(rx1, ox1);
                    ry1 = fast_min(ry1, oy1);
                    rx2 = fast_max(rx2, ox2);
                    ry2 = fast_max(ry2, oy2);
                }
            }
            else {
				overlap[inter_count++].second = j; // intersect
            }
        }        
        
        // Greedily add superpixels according to the intersect overlap
        if (inter_count > 0 && inner_count > 0) {
            // inner overlap
            area_inner = (rx2 - rx1 + 1) * (ry2 - ry1 + 1);
            overlap_max = (float)area_inner / areaI;
            
            // sort straddling set 
            for (j = 0; j < inter_count; j++) {
                // regions combination
				k = overlap[j].second;
				mx1 = fast_min(sp_boxes[k][0], rx1);
				my1 = fast_min(sp_boxes[k][1], ry1);
				mx2 = fast_max(sp_boxes[k][2], rx2);
				my2 = fast_max(sp_boxes[k][3], ry2);
                
                // overlap
                ox1 = fast_max(x1, mx1);
                oy1 = fast_max(y1, my1);
                ox2 = fast_min(x2, mx2);
                oy2 = fast_min(y2, my2);
                
                areaM = (mx2 - mx1 + 1) * (my2 - my1 + 1);
                area_inter = (ox2 - ox1 + 1) * (oy2 - oy1 + 1);
				overlap[j].first = (float)area_inter / (areaM + areaI - area_inter);
            }
            // sort overlaps in descending order
			sort(overlap.begin(), overlap.begin() + inter_count, greater<pair<float, int> >());
			overlap_new = overlap[0].first;
            
            // greedy expansion
            if (overlap_new >= overlap_max) {
                // combine regions
				k = overlap[0].second;
				mx1 = fast_min(sp_boxes[k][0], rx1);
				my1 = fast_min(sp_boxes[k][1], ry1);
				mx2 = fast_max(sp_boxes[k][2], rx2);
				my2 = fast_max(sp_boxes[k][3], ry2);
				j = 1;
                
                while (overlap_new >= overlap_max && j < inter_count) {
                    // update box
                    overlap_max = overlap_new;
                    rx1 = mx1; ry1 = my1; rx2 = mx2; ry2 = my2;
                                        
                    // combine regions
					k = overlap[j].second;
					mx1 = fast_min(sp_boxes[k][0], rx1);
					my1 = fast_min(sp_boxes[k][1], ry1);
					mx2 = fast_max(sp_boxes[k][2], rx2);
					my2 = fast_max(sp_boxes[k][3], ry2);
                    
                    // overlap
                    ox1 = fast_max(x1, mx1);
                    oy1 = fast_max(y1, my1);
                    ox2 = fast_min(x2, mx2);
                    oy2 = fast_min(y2, my2);
                    
                    areaM = (mx2 - mx1 + 1) * (my2 - my1 + 1);
                    area_inter = (ox2 - ox1 + 1) * (oy2 - oy1 + 1);
                    overlap_new = (float)area_inter / (areaM + areaI - area_inter);
					j = j + 1;
                }
            }
        }
        
        // update box
		bool flag = true;
        for (j = 0; j < i; j++) {  // remove duplicate
			if (rx1 == init_boxes[j][0] && ry1 == init_boxes[j][1] && rx2 == init_boxes[j][2] && ry2 == init_boxes[j][3]) {
				init_boxes(i) = -1;
				flag = false;
                break;
            }
        }
        if (flag) {
			init_boxes[i][0] = rx1, init_boxes[i][1] = ry1;
			init_boxes[i][2] = rx2, init_boxes[i][3] = ry2;
            nnzbox ++;
        }   
    }
    
    // sort in descending order w.r.t score
	init_boxes.sort();
    
    // stage 2: straddling expansion
    // We use bounding box instead of segment to compute straddling degree
    // as they yield similar performance while bounding box saves computation
	std::srand(131);//srand((unsigned int)time(NULL));
	vecF inters(numS);
	vecI sp_index(numS);
	for (i = 0; i < nnzbox; i++) {
		int counter = 0;
		rx1 = init_boxes[i][0], ry1 = init_boxes[i][1];
		rx2 = init_boxes[i][2], ry2 = init_boxes[i][3];
        
        for (k = 0; k < numS; k++) {            
            // overlap
			ox1 = fast_max(sp_boxes[k][0], rx1);
			ox2 = fast_min(sp_boxes[k][2], rx2);
			if (ox1 > ox2)  continue;
			oy1 = fast_max(sp_boxes[k][1], ry1);
			oy2 = fast_min(sp_boxes[k][3], ry2);
			if (oy1 > oy2)  continue;

			area_inter = (ox2 - ox1 + 1) * (oy2 - oy1 + 1);
			inter = (float)area_inter / areaS[k];
			inters[counter] = inter;
			sp_index[counter++] = k;
        }

        // add segments with straddling degree >= theta
        for (int t = 0; t < numT; t++) {
			bool flag = false;
            for (j = 0; j < counter; j++) {
				if (inters[j] >= thetas[t]) {
                    // combine regions
                    k = sp_index[j];
					rx1 = fast_min(sp_boxes[k][0], rx1);
					ry1 = fast_min(sp_boxes[k][1], ry1);
					rx2 = fast_max(sp_boxes[k][2], rx2);
					ry2 = fast_max(sp_boxes[k][3], ry2);

                    inters[j] = 0;
					flag = true;
                }
            }

            // final result
            if (flag) {
				Vec4i box(rx1, ry1, rx2, ry2);
				double val = -i * (rand() / double(RAND_MAX));
				out_boxes.pushBack(val, box);
            }
        }
    }
    
    // non maximal suppression, return no more than 5000 proposals
    boxesNms(out_boxes, beta, 5000);
    
    //	release memory
   	delete[] areaS;
}

void SegmentRecursiveBox(vector<Vec4i> &sp_boxes, ValStructVec<float, Vec4i> &init_boxes, ValStructVec<float, Vec4i> &out_boxes, vecF thetas, const float beta, const bool combine)
{
	const int numI = init_boxes.size();
	sort(thetas.rbegin(), thetas.rend());
    
	if (combine)
		out_boxes.reserve(3000);
	else
		out_boxes.reserve(1000);
	SegmentRecursiveBox(sp_boxes, init_boxes, thetas, beta, out_boxes);
	if (combine) {
		for (int i = 0; i < numI; i++)
			out_boxes.pushBack(init_boxes(i), init_boxes[i]);
	}
}
static int **segIdx_imgs;

void initSegMask(int threadNum) 
{
	segIdx_imgs = new int*[threadNum];

	for (int i = 0; i < threadNum; i++) {
		segIdx_imgs[i] = new int[360 * 400];
	}

	int gridr = 4,gridc = 4;
	int rows = 360, cols = 400, segR = ceil(360 / gridr), segC = ceil(400/gridc);
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			int idx = i*cols + j;
			int segi = floor(i / gridr), segj = floor(j / gridc);
			for (int t = 0; t < threadNum; t++) {
				segIdx_imgs[t][idx] = segi * segC + segj;
				segIdx_imgs[t][idx] = segi * segC + segj;
				segIdx_imgs[t][idx] = segi * segC + segj;
				segIdx_imgs[t][idx] = segi * segC + segj;
			}
		}
	}
}

void releaseSegMask(int threadNum)
{
	for (int i = 0; i < threadNum; i++) {
		delete segIdx_imgs[i];
	}
	delete[] segIdx_imgs;
}

inline float diff(Vec3f in1, Vec3f in2)
{
	return sqrt(square(in1[0] - in2[0]) + square(in1[1] - in2[1]) + square(in1[2] - in2[2]));
}

void Segmentation(vector<Vec4i> &boxes, Mat &matIm, float c, int min_size)
{
	const int _w = matIm.cols;
	const int _h = matIm.rows;
	const int _s = _w*_h;
	
	int *idx_imgTemp = segIdx_imgs[omp_get_thread_num()];

	// normalize the indexes of the image in idx_img
	int _max = (_w / 4)*(_h / 4), NUM_IND = 0;
	vecI indexes(_max, 0);
	for (int i = 0; i < _s; i++)
		indexes[idx_imgTemp[i]]++;
	for (int i = 0; i < _max; i++) {
		if (indexes[i] != 0)
			indexes[i] = NUM_IND++;
	}
	int *idx_img = new int[_s];
	for (int i = 0; i < _s; i++)
		idx_img[i] = indexes[idx_imgTemp[i]];
    
	// get adjacent domains
	vector<vecI> adjacent(NUM_IND);
	for (int i = 0; i < NUM_IND; i++)
		adjacent[i].reserve(6);
	int counter = 0;
	const int shift[4] = { 1, _w, -1, -_w };
	for (int row = 1; row < _h - 1; row++) {
		counter = row*_w + 1;
		for (int col = 1; col < _w - 1; col++) {
			int curr = idx_img[counter];
			for (int i = 0; i < 4; i++) {
				int pre = idx_img[counter+shift[i]];
				if (curr > pre&&find(adjacent[curr].begin(), adjacent[curr].end(), pre) == adjacent[curr].end())
					adjacent[curr].push_back(pre);
			}
			counter++;
		}
	}

	// compute the average pixel values of every domain
	counter = 0;
	vector<Vec3f> avgImg(NUM_IND, Vec3f(0, 0, 0));
	vecI pixNum(NUM_IND, 0);
	for (int row = 0; row < _h; row++) {
		Vec3b *ptrIm = matIm.ptr<Vec3b>(row);
		for (int col = 0; col < _w; col++) {
			int value = idx_img[counter++];
			avgImg[value] += ptrIm[col];
			pixNum[value]++;
		}
	}
	int num = 0;
	for (int i = 0; i < NUM_IND; i++) {
		avgImg[i] /= pixNum[i];
		num += adjacent[i].size();
	}

	// compute the edges of adjacent domains
	edge *edges = new edge[num];
	num = 0;
	for (int i = 0; i < NUM_IND; i++) {
		int adjaNum = adjacent[i].size();
		for (int j = 0; j < adjaNum; j++) {
			edges[num].a = i;
			edges[num].b = adjacent[i][j];
			edges[num].w = diff(avgImg[i], avgImg[adjacent[i][j]]);
			num++;
		}
	}

	// segment
	universe *u = segment_graph(NUM_IND, num, edges, c);

	// post process small components
	for (int i = 0; i < num; i++) {
		int a = u->find(edges[i].a);
		int b = u->find(edges[i].b);
		if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
			u->join(a, b);
	}
	int num_css = u->num_sets();
	// compute the indexes of superpixels in idx_img
	num = 1;
	std::memset(indexes.data(), 0, NUM_IND*sizeof(int));
	for (int i = 0; i < NUM_IND; i++) {
		int comp = u->find(i);
		if (!indexes[comp])
			indexes[comp] = num++;
		pixNum[i] = indexes[comp] - 1;
	}
	CV_Assert(num_css == num - 1);

	// get the bounding boxes of the superpixels
	counter = 0;
	boxes.resize(num_css, Vec4i(_w, _h, 1, 1));
	for (int row = 1; row <= _h; row++) {
		for (int col = 1; col <= _w; col++) {
			int curr = pixNum[idx_img[counter++]];
			if (boxes[curr][0] > col)
				boxes[curr][0] = col;
			if (boxes[curr][1] > row)
				boxes[curr][1] = row;
			if (boxes[curr][2] < col)
				boxes[curr][2] = col;
			if (boxes[curr][3] < row)
				boxes[curr][3] = row;
		}
	}

	u->~universe();
	delete []edges;
	delete []idx_img;
}
