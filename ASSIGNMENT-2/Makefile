CXX = g++
CXXFLAGS = -Wall -std=c++11

OBJS = Vector.o Polygon.o Voronoi.o SVG.o main.o
TARGET = voronoi

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJS) $(TARGET)
