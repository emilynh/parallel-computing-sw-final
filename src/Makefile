C = nvcc
NVCCFLAGS = -arch=sm_60 
CFLAGS = -std=c++11

# Thêm tùy chọn -fopenmp vào trình biên dịch C++ của bạn
CFLAGS += -fopenmp

all: align

align: align.cu  
	$(C) $(NVCCFLAGS) $(CFLAGS) -o align.exe align.cu 

clean:
	rm -f align.exe *.dat
