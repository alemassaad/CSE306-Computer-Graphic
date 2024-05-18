CXX = g++
CXXFLAGS = -std=c++11 -O3

OBJS = main_surfaces.o main_light.o vector.o ray.o sphere.o scene.o random.o random_cos.o stb_image_write_impl.o stb_image_impl.o triangle_mesh.o triangle_intersection.o bvh.o bounding_box.o

all: render_surfaces render_light render_rmi

render_surfaces: main_surfaces.o vector.o ray.o sphere.o scene.o random.o random_cos.o stb_image_write_impl.o stb_image_impl.o
	$(CXX) $(CXXFLAGS) -o render_surfaces main_surfaces.o vector.o ray.o sphere.o scene.o random.o random_cos.o stb_image_write_impl.o stb_image_impl.o

render_light: main_light.o vector.o ray.o sphere.o scene.o random.o random_cos.o stb_image_write_impl.o stb_image_impl.o
	$(CXX) $(CXXFLAGS) -o render_light main_light.o vector.o ray.o sphere.o scene.o random.o random_cos.o stb_image_write_impl.o stb_image_impl.o

render_rmi: main_rmi.o vector.o ray.o sphere.o scene.o random.o random_cos.o stb_image_write_impl.o stb_image_impl.o triangle_mesh.o triangle_intersection.o bvh.o bounding_box.o
	$(CXX) $(CXXFLAGS) -o render_rmi main_rmi.o vector.o ray.o sphere.o scene.o random.o random_cos.o stb_image_write_impl.o stb_image_impl.o triangle_mesh.o triangle_intersection.o bvh.o bounding_box.o

main_surfaces.o: main_surfaces.cpp vector.h ray.h sphere.h scene.h stb_image_write.h
	$(CXX) $(CXXFLAGS) -c main_surfaces.cpp

main_light.o: main_light.cpp vector.h ray.h sphere.h scene.h random.h random_cos.h stb_image_write.h
	$(CXX) $(CXXFLAGS) -c main_light.cpp

main_rmi.o: main_rmi.cpp vector.h ray.h sphere.h scene.h random.h random_cos.h stb_image_write.h triangle_mesh.h
	$(CXX) $(CXXFLAGS) -c main_rmi.cpp

vector.o: vector.cpp vector.h
	$(CXX) $(CXXFLAGS) -c vector.cpp

ray.o: ray.cpp ray.h vector.h
	$(CXX) $(CXXFLAGS) -c ray.cpp

sphere.o: sphere.cpp sphere.h ray.h vector.h
	$(CXX) $(CXXFLAGS) -c sphere.cpp

scene.o: scene.cpp scene.h sphere.h ray.h vector.h random_cos.h
	$(CXX) $(CXXFLAGS) -c scene.cpp

random.o: random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp

random_cos.o: random_cos.cpp random_cos.h random.h
	$(CXX) $(CXXFLAGS) -c random_cos.cpp

triangle_mesh.o: triangle_mesh.cpp triangle_mesh.h triangle_intersection.h
	$(CXX) $(CXXFLAGS) -c triangle_mesh.cpp

triangle_intersection.o: triangle_intersection.cpp triangle_intersection.h
	$(CXX) $(CXXFLAGS) -c triangle_intersection.cpp

bvh.o: bvh.cpp bvh.h bounding_box.h
	$(CXX) $(CXXFLAGS) -c bvh.cpp

bounding_box.o: bounding_box.cpp bounding_box.h
	$(CXX) $(CXXFLAGS) -c bounding_box.cpp

stb_image_write_impl.o: stb_image_write_impl.cpp stb_image_write.h
	$(CXX) $(CXXFLAGS) -c stb_image_write_impl.cpp

stb_image_impl.o: stb_image_impl.cpp stb_image.h
	$(CXX) $(CXXFLAGS) -c stb_image_impl.cpp

clean:
	rm -f *.o render_surfaces render_light render_rmi image_surfaces.png image_light.png image_rmi.png
