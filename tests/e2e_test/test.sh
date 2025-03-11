#!/bin/bash

# Проверяем, переданы ли аргументы (1 — путь к бинарнику, 2 — путь к тестам)
if [ $# -lt 2 ]; then
    echo "❌ Использование: $0 <путь_к_бинарнику> <путь_к_папке_с_тестами>"
    exit 1
fi

binary="$1"        # Путь к бинарнику
data_dir="$2"      # Путь к папке с тестами
epsilon=0.001      # Порог сравнения

if [ ! -d "$data_dir" ]; then
    echo "❌ Папка $data_dir не найдена!"
    exit 1
fi

files=("$data_dir"/*.dat)

if [ ${#files[@]} -eq 0 ] || [ "${files[0]}" == "$data_dir/*.dat" ]; then
    echo "❌ Нет файлов с расширением .dat в папке $data_dir."
    exit 1
fi

for file in "${files[@]}"; do
    base_name=$(basename "$file" .dat)        # Получаем имя файла без расширения
    expected_file="$data_dir/$base_name.ans"  # Полный путь к файлу .ans

    echo "📂 Обработка файла: $file"

    if [ ! -f "$expected_file" ]; then
        echo "❌ Файл с ожидаемым результатом ($expected_file) не найден!"
        exit 1
    fi

    echo "🔢 Определитель матрицы:"
    result=$("$binary" < "$file")     # Запускаем программу и сохраняем результат
    expected=$(cat "$expected_file")  # Читаем ожидаемый результат

    echo "🔹 Результат: $result"
    echo "🔹 Ожидалось: $expected"

    # Разделитель дробной части — точка
    export LC_NUMERIC="en_US.UTF-8"

    # Вычисление разницы и сравнение с epsilon
    difference=$(awk -v r="$result" -v e="$expected" 'BEGIN { print ((r - e) < 0) ? -(r - e) : (r - e) }')

    if [ "$(echo "$difference < $epsilon" | bc -l)" -eq 1 ]; then
        echo "✅ Тест пройден!"
    else
        echo "❌ Ошибка! Разница $difference превышает порог $epsilon."
        exit 1
    fi
    echo "-----------------------------------"
done
