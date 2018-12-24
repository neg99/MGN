# MGN
Modified Gauss-Newton algorithm

# Код на R, использующийся в кандидатской

**Прежде, чем запускать код, необходимо установить вспомогательный пакет из каталога [rhlra2](rhlra2/)**.

Необходимо установить следующие пакеты:
`install.packages(c("fftw", "svd", "Matrix", "orthopolynom"))`

## Описание каталогов

### [auxiliary](auxiliary/)

Содержит код методов, использующийся в моделировании (Cadzow, MGN, VP, поиск весов)

### [plots](plots/)

Код для построения сравнения методов вычисления базисов, а так же сравнения методов
MGN/VPGN. Результат - рисунки (в кандидатской/статье).

## Описание файлов

Общее правило - файл достаточно просто открыть, установить в RStudio рабочий каталог в место файла,
и запустить его целиком. Дальше он сам сделает `source` для всех зависимостей.

Для того, чтобы сразу начать, можно открыть и почитать смысл параметров методов в [ридми в auxiliary](auxiliary).
А прям в этом каталоге есть несколько примеров решения задачи HSLRA для белого/красного шума. Их можно сразу открыть и запустить.