CXX = g++
CXXFLAGS = -Wall -Werror -Wextra -std=c++17
LDFLAGS = -lgtest -lgtest_main -pthread
TARGET = test_add
LIBRARY = libS21Matrix.a
SRCS = s21_matrix_oop.cpp test.cpp
OBJS = $(SRCS:.cpp=.o)
all: $(LIBRARY) $(TARGET)

# Цель для создания статической библиотеки
$(LIBRARY): s21_matrix_oop.o
	ar rcs $@ $^

# Цель для компиляции тестов
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Правило для компиляции объектных файлов
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Если бинарник в текущей дирректории, не компилируем еще раз
test: $(TARGET)
	@if [ -f $(TARGET) ]; then \
		./$(TARGET) || true; \
	fi

memory: $(TARGET)
	@if [ -f $(TARGET) ]; then \
		valgrind --leak-check=full --track-origins=yes ./$(TARGET); \
	fi

clean:
	rm -rf test_add .vscode cover* *.gcda *.gcno s21_matrix_oop.h.gch *.o *.a

.PHONY: all clean

# gcov_report: clean
# 	g++ -fprofile-arcs -ftest-coverage -o test_add s21_* test.cpp -lgtest -lgtest_main -pthread
# 	./test_add
# 	gcov *.cpp
# 	geninfo -o coverage.info .
# 	genhtml coverage.info -o coverage_report

gcov_report: clean
	$(CXX) -fprofile-arcs -ftest-coverage $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)
	./$(TARGET) # Запускаем тесты
	gcov s21_matrix_oop.cpp 
	gcov s21_matrix_oop.h   
	geninfo -o coverage.info .
	genhtml coverage.info -o coverage_report


clang:
	clang-format --style=Google -n *.cpp *.h

# docker targets
OS := $(shell uname -s)
install_docker:
	@if [ "$(OS)" = "Linux" ]; then \
		sudo apt-get update; \
		sudo apt-get install -y apt-transport-https ca-certificates curl software-properties-common; \
		curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -; \
		sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $$(lsb_release -cs) stable"; \
		sudo apt-get update; \
		sudo apt-get install -y docker-ce; \
	elif [ "$(OS)" = "Darwin" ]; then \
		brew install --cask docker; \

build: stop
	docker build -t s21_matrix .
	docker run -d -it --name s21_matrix s21_matrix /bin/sh

run:
	docker exec -it s21_matrix /bin/sh

stop:
	docker stop s21_matrix || true
	docker rm s21_matrix --force || true 
	docker image rm s21_matrix --force || true 

restart: stop build

refresh:
	docker exec s21_matrix ls /app
	docker cp . s21_matrix:/app

pull_report:
	docker cp s21_matrix:/app/coverage_report ./coverage_report
