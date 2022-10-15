1. Сборка
=========

1.1 Windows, VisualStudio
-------------------------
Проверено на Windows 10, Visual Studio 2019

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

- cоздать папку build в корне проекта (папка ns/ репозитория)

- скопировать скрипт winbuild64.bat в папку build. Далее вносить изменения
  только в скопированном файле.

- скрипт написан для версии Visual Studio 2019. Если используется другая версия,
  изменить в скрипте значения переменных CMGenerator и CMToolset на соответствующие вашей версии.

- в скрипте winbuild64.bat подставить пути к распакованным папкам в соответствующие строчки.
  Например

SET CMBOOST_ROOT="c:/lib/boost_1_78_0"
SET CMAMGCL_ROOT="c:/lib/amgcl-1.4.2"

- запустить скрипт winbuild64.bat из папки build

- После сборки в папке build появится проект VisualStudio.
  Сборка бинарных файлов будет осуществлятся в папку build/bin.
  Выходные файлы при использовании функции from_output_path будут писаться в build/output

1.2 Linux, GCC
--------------

- установить boost и amgcl в систему
- в консоли
mkdir build
cd build
cmake ..
make
- программы соберутся в папку build/bin
