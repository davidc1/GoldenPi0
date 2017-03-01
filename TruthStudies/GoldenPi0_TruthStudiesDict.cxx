// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GoldenPi0_TruthStudiesDict

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
#include "BNBPi0Properties.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLBNBPi0Properties_Dictionary();
   static void larlitecLcLBNBPi0Properties_TClassManip(TClass*);
   static void *new_larlitecLcLBNBPi0Properties(void *p = 0);
   static void *newArray_larlitecLcLBNBPi0Properties(Long_t size, void *p);
   static void delete_larlitecLcLBNBPi0Properties(void *p);
   static void deleteArray_larlitecLcLBNBPi0Properties(void *p);
   static void destruct_larlitecLcLBNBPi0Properties(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::BNBPi0Properties*)
   {
      ::larlite::BNBPi0Properties *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::BNBPi0Properties));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::BNBPi0Properties", "BNBPi0Properties.h", 25,
                  typeid(::larlite::BNBPi0Properties), DefineBehavior(ptr, ptr),
                  &larlitecLcLBNBPi0Properties_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::BNBPi0Properties) );
      instance.SetNew(&new_larlitecLcLBNBPi0Properties);
      instance.SetNewArray(&newArray_larlitecLcLBNBPi0Properties);
      instance.SetDelete(&delete_larlitecLcLBNBPi0Properties);
      instance.SetDeleteArray(&deleteArray_larlitecLcLBNBPi0Properties);
      instance.SetDestructor(&destruct_larlitecLcLBNBPi0Properties);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::BNBPi0Properties*)
   {
      return GenerateInitInstanceLocal((::larlite::BNBPi0Properties*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::BNBPi0Properties*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLBNBPi0Properties_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::BNBPi0Properties*)0x0)->GetClass();
      larlitecLcLBNBPi0Properties_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLBNBPi0Properties_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLBNBPi0Properties(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::BNBPi0Properties : new ::larlite::BNBPi0Properties;
   }
   static void *newArray_larlitecLcLBNBPi0Properties(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::BNBPi0Properties[nElements] : new ::larlite::BNBPi0Properties[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLBNBPi0Properties(void *p) {
      delete ((::larlite::BNBPi0Properties*)p);
   }
   static void deleteArray_larlitecLcLBNBPi0Properties(void *p) {
      delete [] ((::larlite::BNBPi0Properties*)p);
   }
   static void destruct_larlitecLcLBNBPi0Properties(void *p) {
      typedef ::larlite::BNBPi0Properties current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::BNBPi0Properties

namespace {
  void TriggerDictionaryInitialization_libGoldenPi0_TruthStudies_Impl() {
    static const char* headers[] = {
"BNBPi0Properties.h",
0
    };
    static const char* includePaths[] = {
"/uboone/app/users/davidc1/LArLite/core",
"/grid/fermiapp/products/larsoft/root/v6_04_02/Linux64bit+2.6-2.12-e7-prof/include",
"/uboone/app/users/davidc1/LArLite/UserDev/GoldenPi0/TruthStudies/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$BNBPi0Properties.h")))  BNBPi0Properties;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "BNBPi0Properties.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::BNBPi0Properties", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGoldenPi0_TruthStudies",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGoldenPi0_TruthStudies_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGoldenPi0_TruthStudies_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGoldenPi0_TruthStudies() {
  TriggerDictionaryInitialization_libGoldenPi0_TruthStudies_Impl();
}
