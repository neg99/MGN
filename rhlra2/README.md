# Вспомогательный пакет rhlra2

## Инструкция для linux/macOS

Требования:
1. Окружение для сборки (работающий компилятор `GCC`/`clang` и `dev`-пакеты для `R`, если такие есть).
2. `fftw` (устанавливается через пакетный менеджер (`apt-get`, `brew`, etc.)).

Установка (несколько вариантов):
1. Через RStudio стандартным способом (New project -> Existing directory -> указать директорию с rhlra. Затем кнопка `Clean and rebuild`(или `Build and reload`) для сборки. Лежит в правом верхнем углу окна во вкладке `Build`).
2. Из консоли: в каталоге `rhlra` выполнить `R CMD INSTALL .`

## Инструкция для Windows

Требования:
1. [Rtools](https://cran.r-project.org/bin/windows/Rtools/). Желательно установить в `C:\RBuildTools`. 
2. `fftw` необходимо взять [отсюда](http://www.fftw.org/install/windows.html). Скачать нужную версию собранной библиотеки из раздела "Precompiled FFTW 3.3.5 Windows DLLs", **обязательно** распаковать всё содержимое в `C:\RBuildTools\fftw` и **прописать этот каталог в PATH** (см. инструкцию для Windows 10, например, [здесь](https://superuser.com/questions/949560/how-do-i-set-system-environment-variables-in-windows-10); **обязательно** после редактирования PATH необходимо делать логаут/логин (или перезагрузить систему, _ох уж этот шиндовс..._)).
3. Пакет `devtools` для `R` (`install.packages("devtools")`).

Установка (пошагово):
1. Открыть RStudio, в ней загрузить пакет `devtools` (`library(devtools)`).
2. Выполнить `find_rtools()`, убедиться, что она возвращает `TRUE`.
3. Далее через RStudio стандартным способом (New project -> Existing directory -> указать директорию с rhlra. Затем кнопка `Clean and rebuild`(или `Build and reload`) для сборки. Лежит в правом верхнем углу окна во вкладке `Build`).