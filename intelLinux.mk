# Version 	: 0.0.1
# Date		: 26.01.2012
# Author 	: Cedric.Bilat@he-arc.ch
#
# Contrainte
#
#	Dans les definitions de variables ci-dessous, ne laisser jamais d'espace a la fin!!!
#

#########
#public #
#########

SRC_AUX+=../PRODUCTION/commonExt/include
SRC_AUX+=../PRODUCTION/commonExt/linux

#Option de compilation 
#Elles seront ensuite surchargï¿½es automatiquement (par exemple en fonction du type de la target)
#Compilation: (minimal conseillé : vide) 
#CXXFLAGS:= -w 
CXXFLAGS+= -axSSSE3
CXXFLAGS+= -fast -O3 -parallel -par-threshold100  
CXXFLAGS+= -openmp

#Option de compilation Static Lib 
ARFLAGS:=

#Option de linkage
#Elles seront ensuite surchargï¿½es automatiquement (par exemple en fonction du type de la target)
#Version minimale : Ne rien metter laisser vierge!
LDFLAGS_AUX:=

#Injection de variable dans le code (same as #define XXX YYY)
CODE_DEFINE_VARIABLES:=

#separateur est espace
EXCLUDE_LIBRARY_FILES:=
ADD_LIBRARY_FILES:=

#########
#Private#
#########

COMPILATEUR:=INTEL
OS:=Linux

ROOT_MAKEFILE_PRIVATE:=../PRODUCTION/MANAGER/makefile/private
-include dataProject.mk
-include $(ROOT_MAKEFILE_PRIVATE)/fileSet.mk
-include $(ROOT_MAKEFILE_PRIVATE)/flags_GccIntel.mk
-include $(ROOT_MAKEFILE_PRIVATE)/makeFileCPP.mk
 
 
#########
#  End  #
#########