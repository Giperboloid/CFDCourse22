1. Сборка
=========

1.1 Windows, VisualStudio
-------------------------

- Cкачать amgcl
https://github.com/ddemidov/amgcl/archive/refs/tags/1.4.2.zip

- Скачать boost 
https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.zip

- установить cmake
https://github.com/Kitware/CMake/releases/download/v3.24.2/cmake-3.24.2-windows-x86_64.msi

- Распаковать эти архивы. Например в C:/lib/.
Чтобы в этой папке были директории c содержанием архивов
amgcl-1.4.2\*      {внутри неё должна быть в том числе папка amgcl\} 
boost_1_78_0\*     {внутри неё должна быть в том числе папка boost\}

- в файле winbuild64.bat подставить пути к этим папкам в соответствующие строчки. Например

SET CMBOOST_ROOT="c:/lib/boost_1_78_0"
SET CMBAMGCL_ROOT="c:/lib/amgcl-1.4.2"

- запустить скрипт winbuild64.bat. Появится папка build64 с проектом. 

- После сборки в папке build64/bin появятся исполняемые файлы

1.2 Linux, GCC
--------------

- установить boost и amgcl в систему
- в консоли
mkdir build
cd build
cmake ..
make
- программы соберутся в папку build/bin
