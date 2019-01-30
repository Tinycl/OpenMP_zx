#include <windows.h>
#include <iostream>
#include <omp.h>

using namespace std;

void SimpleWithOMP(int n, float *a, float *b)
{
	int i;
	int iam, nthd, iend, istart;
	#pragma omp parallel for //default(shared) private(iam,nthd, istart, iend)
	for (i = 1; i < n; i++)
	{
		b[i] = (a[i] + a[i - 1]) / 2.0 * 3.0;
		//printf("threadid is %d  i is %d \n",omp_get_thread_num(),i); //打印core number
	}
}
void SimpleNoOMP( int n, float *a, float *b)
{
	int i;
	for (i = 1; i < n; i++)
	{
		b[i] = (a[i] + a[i - 1]) / 2.0 * 3.0;
	}
}

void subdomain(float *x, int istart, int ipoints)
{
	int i;
	for (i = 0; i < ipoints; i++)
	{
		x[istart + i] = 123.456;
		printf("threadid is %d \n",omp_get_thread_num());
	}
}
//#pragma omp parallel 实现 for 块的并行
void sub(float *x, int npoints)
{
	int iam, nt, ipoints, istart;
	#pragma omp parallel default(shared) private(iam,nt, ipoints, istart)
	{
		iam = omp_get_thread_num(); //哪个core执行  可以认为core数和可以并行优化的线程数一致
		nt = omp_get_num_threads(); //总的core个数
		ipoints = npoints / nt; //size of partition 拆分数组
		istart = iam * ipoints; //starting arrray index
		if (iam == nt - 1)   //last thread may do more or less
		{
			ipoints = npoints - istart;
		}
		subdomain(x,istart,ipoints);
	}
}
//parallel 实现 数组 for 遍历计算 等价于 parallel for
void array_iterate_cal(float *array, int length)
{
	int iam, istart, iend, nthrd, iter;
	#pragma omp parallel default(shared) private(iam,istart, iend, nthrd, iter) num_threads(4)
	{
		iam = omp_get_thread_num();  //并行区域 正在执行的线程
		nthrd = omp_get_num_threads(); //并行区域 共多少个线程
		iend = length / nthrd;   //把数组元素 分给每个线程处理，每一段的长度
		istart = iam * iend;    //每一段的开始
		if (iam == nthrd - 1)  //最后一个线程，处理的元素可能少 小于iend
		{
			iend = length - istart;
		}
		for (iter = 0; iter < iend; iter++)
		{
			array[istart + iter] = 1.0*2.0 + 1.0;
			printf("threadid is %d  iter is %d \n", omp_get_thread_num(), istart + iter);
		}
	}
}

void work1()
{

}
void work2()
{

}
static omp_lock_t g_lock;
int main()
{
	int arraysize = 10;// 2592 * 1944;// 2592 * 1944; //2048*1536
	float *a = new float[arraysize];
	for (int i = 0; i < arraysize; i++)
	{
		a[0] = i * 1.0;
	}
	float *b = new float[arraysize];
#ifndef _OPENMP
	fprintf(stderr, "OpenMP not supported");
#endif
	omp_set_num_threads(4); //设置4个core并行处理 4core cpu
	//测试优化性能
	LARGE_INTEGER timestart;
	LARGE_INTEGER timeend;
	LARGE_INTEGER frequency;
	double quadpart = 0;
	double elapsed = 0;
	printf("\/************optimation difference************\/\n");
	
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&timestart);
	quadpart = (double)frequency.QuadPart; 
	SimpleWithOMP(arraysize,a,b);
	QueryPerformanceCounter(&timeend);
	elapsed = (timeend.QuadPart - timestart.QuadPart) / quadpart; 
	printf("omp optimation for loop NO : %f s \n", elapsed);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&timestart);
	quadpart = (double)frequency.QuadPart;
	SimpleNoOMP(arraysize, a, b);
	QueryPerformanceCounter(&timeend);
	elapsed = (timeend.QuadPart - timestart.QuadPart) / quadpart;
	printf("omp optimation for loop YES : %f s \n", elapsed);
	if (a != NULL) {
		delete[] a;
	}
	if (b != NULL) {
		delete[] b;
	}
	
	//任务分配，数组拆分
	printf("\/************array partition************\/\n");
	float array[10];
	sub(array, 10);
	//嵌套
	printf("\/************nest test************\/\n");
	float *array1 = new float[2555*2555];
	double ompstarttime = 0, ompendtime = 0;
	omp_set_nested(1); //设置一层并行嵌套
	ompstarttime = omp_get_wtime();
	#pragma omp parallel for
	for (int i = 0; i < 2555; i++) {
		//#pragma omp parallel for 嵌套一层循环时， 只需在最外层加指令，否则结果负优化（串行）
		for(int j = 0; j < 2555; j++){
			array1[i*2555+j] = 2.0*1.5 + 1;
		}
	}
	ompendtime = omp_get_wtime();
	printf("use 1 nest to cal 2 dim : %f s \n", ompendtime - ompstarttime);
	
	omp_set_nested(0); //不设置并行嵌套
	ompstarttime = omp_get_wtime();
	for (int i = 0; i < 2555; i++) {
		for (int j = 0; j < 2555; j++) {
			array1[i*2555+j] = 2.0*1.5 + 1;
		}
	}
	ompendtime = omp_get_wtime();
	printf("don't use nest to cal 2 dim : %f s \n", ompendtime - ompstarttime);
	delete[] array1;
	//设置动态调整
	printf("\/************set dynamic************\/\n");
	float array2[10];
	omp_set_dynamic(0); //禁用动态调整
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < 10; i++) {
		array2[i] = 1.0*2.0 + 1.0;
		printf("threadid is %d \n", omp_get_thread_num());
	}

	//线程同步之omp_lock
	printf("\/************mutex lock************\/\n");
	omp_init_lock(&g_lock);  //初始化lock
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < 10; i++) {
		omp_set_lock(&g_lock);  //获得互斥锁
		printf("threadid is %d +\n", omp_get_thread_num()); //下面3个printf不会打断
		printf("threadid is %d -\n", omp_get_thread_num());
		printf("threadid is %d *\n", omp_get_thread_num());
		omp_unset_lock(&g_lock);  //释放互斥所
	}
	omp_destroy_lock(&g_lock);

	//nowait
	printf("\/************nowait************\/\n");
	#pragma omp parallel 
	{
		//pragma omp for 默认后面有个隐式Barrier
		#pragma omp for nowait     // 两个for block可以并行执行即 +-+-- 而不是等到第一个for完成后执行下一个for +++---
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d +\n", omp_get_thread_num()); 
		}
		#pragma omp for 
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d -\n", omp_get_thread_num());
		}
	}
	
	//barrier
	printf("\/************barrier************\/\n");
	#pragma omp parallel 
	{
		
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d +\n", omp_get_thread_num());
		}
		#pragma omp barrier  //不显示调用barrier， 两个for块并行执行。显示调用barrier，先执行完第一个for，在执行第二个for
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d -\n", omp_get_thread_num());
		}
	}

	//master
	printf("\/************master************\/\n");
	#pragma omp parallel 
	{
		#pragma omp master  //该for并行块只有主线程执行
		{
			for (int i = 0; i < 10; i++) {
				printf("threadid is %d +\n", omp_get_thread_num()); //打印10次
			}
		}
		
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d -\n", omp_get_thread_num()); //打印10*core次
		}
	}

	//pareallel -> parallel for
	printf("\/************parallel -> parallel for************\/\n");
	array_iterate_cal(array2,10);

	//sections
	printf("\/************sections************\/\n");
	#pragma omp parallel sections
	{
		//两个section for 是并行执行的
		#pragma omp section  //只有一个线程printf10次 
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d +\n", omp_get_thread_num());
		}
		#pragma omp section //只有另外一个线程printf10次
		for (int i = 0; i < 10; i++) {
			printf("threadid is %d -\n", omp_get_thread_num());
		}
	}
	//atomic
	printf("\/************atomic************\/\n");
	int sum = 0;
	#pragma omp parallel for
	for (int i = 1; i <= 50; i++)
	{
		//#pragma omp atomic
		sum += i;
		printf("threadid is %d \n", omp_get_thread_num());
	}
	
	//single
	printf("\/************single************\/\n");
	#pragma omp parallel
	{
		#pragma omp single  //并行代码域内，只有一个线程，串行执行该代码
		printf("Beginning work1 \n");
		work1();
		#pragma omp single
		printf("Finishing work1 \n");

		#pragma omp single nowait
		printf("Finish work1 and beginning work2 \n");
		work2();
	}
	system("pause");
	return 0;
}