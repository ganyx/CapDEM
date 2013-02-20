CONFIG += warn_on debug -fopenmp
DEPENDPATH += . src src/post_process
HEADERS += src/cell.h \
 src/config.h \
 src/contact.h \
 src/inout.h \
 src/matrix.h \
 src/mesh.h \
 src/parameter.h \
 src/particle.h \
 src/profile.h \
 src/run.h \
 src/soft_dynamics.h \
 src/stat.h \
 src/post_process.h
INCLUDEPATH += . src src/post_process
LIBS += -fopenmp
SOURCES += src/cdinit.cpp \
 src/cell.cpp \
 src/config.cpp \
 src/contact.cpp \
 src/inout.cpp \
 src/main.cpp \
 src/matrix.cpp \
 src/mesh.cpp \
 src/parameter.cpp \
 src/particle.cpp \
 src/profile.cpp \
 src/run.cpp \
 src/stat.cpp \
 src/post_process.cpp
TARGET =
TEMPLATE = app
