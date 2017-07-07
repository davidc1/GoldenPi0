// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GoldenPi0_Pi0SimchStudiesDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "Pi0HitThresholdStudies.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLPi0HitThresholdStudies_Dictionary();
   static void larlitecLcLPi0HitThresholdStudies_TClassManip(TClass*);
   static void *new_larlitecLcLPi0HitThresholdStudies(void *p = 0);
   static void *newArray_larlitecLcLPi0HitThresholdStudies(Long_t size, void *p);
   static void delete_larlitecLcLPi0HitThresholdStudies(void *p);
   static void deleteArray_larlitecLcLPi0HitThresholdStudies(void *p);
   static void destruct_larlitecLcLPi0HitThresholdStudies(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::Pi0HitThresholdStudies*)
   {
      ::larlite::Pi0HitThresholdStudies *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::Pi0HitThresholdStudies));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::Pi0HitThresholdStudies", "Pi0HitThresholdStudies.h", 29,
                  typeid(::larlite::Pi0HitThresholdStudies), DefineBehavior(ptr, ptr),
                  &larlitecLcLPi0HitThresholdStudies_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::Pi0HitThresholdStudies) );
      instance.SetNew(&new_larlitecLcLPi0HitThresholdStudies);
      instance.SetNewArray(&newArray_larlitecLcLPi0HitThresholdStudies);
      instance.SetDelete(&delete_larlitecLcLPi0HitThresholdStudies);
      instance.SetDeleteArray(&deleteArray_larlitecLcLPi0HitThresholdStudies);
      instance.SetDestructor(&destruct_larlitecLcLPi0HitThresholdStudies);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::Pi0HitThresholdStudies*)
   {
      return GenerateInitInstanceLocal((::larlite::Pi0HitThresholdStudies*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::Pi0HitThresholdStudies*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLPi0HitThresholdStudies_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::Pi0HitThresholdStudies*)0x0)->GetClass();
      larlitecLcLPi0HitThresholdStudies_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLPi0HitThresholdStudies_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLPi0HitThresholdStudies(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::Pi0HitThresholdStudies : new ::larlite::Pi0HitThresholdStudies;
   }
   static void *newArray_larlitecLcLPi0HitThresholdStudies(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::Pi0HitThresholdStudies[nElements] : new ::larlite::Pi0HitThresholdStudies[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLPi0HitThresholdStudies(void *p) {
      delete ((::larlite::Pi0HitThresholdStudies*)p);
   }
   static void deleteArray_larlitecLcLPi0HitThresholdStudies(void *p) {
      delete [] ((::larlite::Pi0HitThresholdStudies*)p);
   }
   static void destruct_larlitecLcLPi0HitThresholdStudies(void *p) {
      typedef ::larlite::Pi0HitThresholdStudies current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::Pi0HitThresholdStudies

namespace {
  void TriggerDictionaryInitialization_libGoldenPi0_Pi0SimchStudies_Impl() {
    static const char* headers[] = {
"Pi0HitThresholdStudies.h",
0
    };
    static const char* includePaths[] = {
"/a/share/westside/dcaratelli/larlite/core",
"/a/share/westside/dcaratelli/larlite/UserDev/RecoTool",
"/a/share/westside/dcaratelli/larlite/UserDev/RecoTool/FANN/include",
"/a/apps/local/root-6.04.00/include",
"/a/share/westside/dcaratelli/larlite/UserDev/GoldenPi0/Pi0SimchStudies/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$Pi0HitThresholdStudies.h")))  Pi0HitThresholdStudies;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "Pi0HitThresholdStudies.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::Pi0HitThresholdStudies", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGoldenPi0_Pi0SimchStudies",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGoldenPi0_Pi0SimchStudies_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGoldenPi0_Pi0SimchStudies_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGoldenPi0_Pi0SimchStudies() {
  TriggerDictionaryInitialization_libGoldenPi0_Pi0SimchStudies_Impl();
}
