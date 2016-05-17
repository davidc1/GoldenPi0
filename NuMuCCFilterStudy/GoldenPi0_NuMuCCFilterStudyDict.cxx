// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GoldenPi0_NuMuCCFilterStudyDict

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
#include "SearchPFPartHierarchy.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLSearchPFPartHierarchy_Dictionary();
   static void larlitecLcLSearchPFPartHierarchy_TClassManip(TClass*);
   static void *new_larlitecLcLSearchPFPartHierarchy(void *p = 0);
   static void *newArray_larlitecLcLSearchPFPartHierarchy(Long_t size, void *p);
   static void delete_larlitecLcLSearchPFPartHierarchy(void *p);
   static void deleteArray_larlitecLcLSearchPFPartHierarchy(void *p);
   static void destruct_larlitecLcLSearchPFPartHierarchy(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::SearchPFPartHierarchy*)
   {
      ::larlite::SearchPFPartHierarchy *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::SearchPFPartHierarchy));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::SearchPFPartHierarchy", "SearchPFPartHierarchy.h", 25,
                  typeid(::larlite::SearchPFPartHierarchy), DefineBehavior(ptr, ptr),
                  &larlitecLcLSearchPFPartHierarchy_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::SearchPFPartHierarchy) );
      instance.SetNew(&new_larlitecLcLSearchPFPartHierarchy);
      instance.SetNewArray(&newArray_larlitecLcLSearchPFPartHierarchy);
      instance.SetDelete(&delete_larlitecLcLSearchPFPartHierarchy);
      instance.SetDeleteArray(&deleteArray_larlitecLcLSearchPFPartHierarchy);
      instance.SetDestructor(&destruct_larlitecLcLSearchPFPartHierarchy);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::SearchPFPartHierarchy*)
   {
      return GenerateInitInstanceLocal((::larlite::SearchPFPartHierarchy*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::SearchPFPartHierarchy*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLSearchPFPartHierarchy_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::SearchPFPartHierarchy*)0x0)->GetClass();
      larlitecLcLSearchPFPartHierarchy_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLSearchPFPartHierarchy_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLSearchPFPartHierarchy(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::SearchPFPartHierarchy : new ::larlite::SearchPFPartHierarchy;
   }
   static void *newArray_larlitecLcLSearchPFPartHierarchy(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::SearchPFPartHierarchy[nElements] : new ::larlite::SearchPFPartHierarchy[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLSearchPFPartHierarchy(void *p) {
      delete ((::larlite::SearchPFPartHierarchy*)p);
   }
   static void deleteArray_larlitecLcLSearchPFPartHierarchy(void *p) {
      delete [] ((::larlite::SearchPFPartHierarchy*)p);
   }
   static void destruct_larlitecLcLSearchPFPartHierarchy(void *p) {
      typedef ::larlite::SearchPFPartHierarchy current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::SearchPFPartHierarchy

namespace {
  void TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy_Impl() {
    static const char* headers[] = {
"SearchPFPartHierarchy.h",
0
    };
    static const char* includePaths[] = {
"/home/david/uboone/larlite/core",
"/home/david/SOFTWARE/ROOT_v60410/include",
"/home/david/uboone/larlite/UserDev/GoldenPi0/NuMuCCFilterStudy/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$SearchPFPartHierarchy.h")))  SearchPFPartHierarchy;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "SearchPFPartHierarchy.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::SearchPFPartHierarchy", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGoldenPi0_NuMuCCFilterStudy",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy() {
  TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy_Impl();
}
