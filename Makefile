NAME := test
#TARGET := readcorsika
TARGET := readtrack

KM2AVer := 
HOST :=$(shell hostname)

OBJDIR :=.
OBJS := $(OBJDIR)/main.o

DEFINES  := -I. -I$(OBJDIR) `root-config --cflags`
LDFLAGS  := `root-config --libs`
LDFLAGS  += -lvdt

CXXFLAGS := -O3 -fPIC

$(NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/main.o: $(TARGET).C
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

debug: CXXFLAGS := $(CXXFLAGS) -g
debug: $(NAME)

clean:
	rm -fv $(OBJDIR)/main.o $(NAME)
