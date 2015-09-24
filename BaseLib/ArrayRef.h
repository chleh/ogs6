#pragma once

#include <cstddef>
#include <stdexcept>

namespace BaseLib
{

template<typename T, std::size_t... Sizes>
class ArrayRef;


/**
 * Dynamically sized non-owning C-array wrapper
 */
template<typename T>
class ArrayRef<T>
{
public:
    // Type
    using value_type = T;
    using size_type  = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using iterator = value_type*;
    using const_iterator = const value_type*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    explicit ArrayRef(T* data) : _data{data} {}

    ArrayRef(ArrayRef<T> const&) = default; // all copies of an instance will point to the __same__ data!
    ArrayRef(ArrayRef<T>&&) = default;      // move will actually copy this reference

    /* operator= might be confusing if used in actual code */
    ArrayRef<T>& operator=(ArrayRef<T> const&) = default;
    ArrayRef<T>& operator=(ArrayRef<T> &&) = default;

    iterator        begin()       { return _data; }
    const_iterator  begin() const { return _data; }
    const_iterator cbegin() const { return _data; }
    iterator          end()       { return _data + _size; }
    const_iterator    end() const { return _data + _size; }
    const_iterator   cend() const { return _data + _size; }

    reverse_iterator        rbegin()       { return _data + _size - 1; }
    const_reverse_iterator  rbegin() const { return _data + _size - 1; }
    const_reverse_iterator crbegin() const { return _data + _size - 1; }
    reverse_iterator          rend()       { return _data - 1; }
    const_reverse_iterator    rend() const { return _data - 1; }
    const_reverse_iterator   crend() const { return _data - 1; }

    const_pointer data() const { return _data; }
    pointer       data()       { return _data; }

    reference                 front()       { return *_data; }
    constexpr const_reference front() const { return *_data; }
    reference                  back()       { return _data[_size-1]; }
    constexpr const_reference  back() const { return _data[_size-1]; }

    reference       operator[](std::size_t pos)       { return _data[pos]; }
    const_reference operator[](std::size_t pos) const { return _data[pos]; }

    reference at(std::size_t pos) {
        if (pos >= _size) throw std::out_of_range("Given index exceeds array size.");
        return _data[pos];
    }
    const_reference at (std::size_t pos) const {
        if (pos >= _size) throw std::out_of_range("Given index exceeds array size.");
        return _data[pos];
    }

    constexpr bool empty() const { return _size == 0; }
    constexpr size_type size() const { return _size; }
    constexpr size_type max_size() const { return _size; }

    void fill(const T& value)     { for (auto& v : *this) v = value; }
    void swap(ArrayRef<T>& other) {
        std::swap(_size, other._size);
        std::swap(_data, other._data);
    }

private:
    size_type _size;
    T* _data;
};



template<typename T, std::size_t Size>
class ArrayRef<T, Size>
{
public:
    // Type
    using value_type = T;
    using size_type  = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using iterator = value_type*;
    using const_iterator = const value_type*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    explicit ArrayRef(T* data) : _data{data} {}

    ArrayRef(ArrayRef<T, Size> const&) = default; // all copies of an instance will point to the __same__ data!
    ArrayRef(ArrayRef<T, Size>&&) = default;      // move will actually copy this reference

    /* operator= might be confusing if used in actual code */
    ArrayRef<T, Size>& operator=(ArrayRef<T, Size> const&) = delete;  // all copies of an instance will point to the __same__ data!
    ArrayRef<T, Size>& operator=(ArrayRef<T, Size> &&) = delete;      // move will actually copy this reference

    iterator        begin()       { return _data; }
    const_iterator  begin() const { return _data; }
    const_iterator cbegin() const { return _data; }
    iterator          end()       { return _data + Size; }
    const_iterator    end() const { return _data + Size; }
    const_iterator   cend() const { return _data + Size; }

    reverse_iterator        rbegin()       { return _data + Size - 1; }
    const_reverse_iterator  rbegin() const { return _data + Size - 1; }
    const_reverse_iterator crbegin() const { return _data + Size - 1; }
    reverse_iterator          rend()       { return _data - 1; }
    const_reverse_iterator    rend() const { return _data - 1; }
    const_reverse_iterator   crend() const { return _data - 1; }

    const_pointer data() const { return _data; }
    pointer       data()       { return _data; }

    reference                 front()       { return *_data; }
    constexpr const_reference front() const { return *_data; }
    reference                  back()       { return _data[Size-1]; }
    constexpr const_reference  back() const { return _data[Size-1]; }

    reference       operator[](std::size_t pos)       { return _data[pos]; }
    const_reference operator[](std::size_t pos) const { return _data[pos]; }

    reference at(std::size_t pos) {
        if (pos >= Size) throw std::out_of_range("Given index exceeds array size.");
        return _data[pos];
    }
    const_reference at (std::size_t pos) const {
        if (pos >= Size) throw std::out_of_range("Given index exceeds array size.");
        return _data[pos];
    }

    constexpr bool empty() const { return Size == 0; }
    constexpr size_type size() const { return Size; }
    constexpr size_type max_size() const { return Size; }

    void fill(const T& value)           { for (auto& v : *this) v = value; }
    void swap(ArrayRef<T, Size>& other) { std::swap(_data, other._data); }

private:
    T* _data;
};

template<typename T, std::size_t Size>
bool operator== (const ArrayRef<T, Size>& lhs, const ArrayRef<T, Size> rhs)
{
    for (auto i=decltype(lhs.size()){0}; i<lhs.size(); ++i) {
        if (lhs[i] != rhs[i]) return false;
    }
    return true;
}

template<typename T, std::size_t Size>
bool operator!= (const ArrayRef<T, Size>& lhs, const ArrayRef<T, Size> rhs)
{
    return ! (lhs == rhs);
}

template<typename T, std::size_t Size>
bool operator<= (const ArrayRef<T, Size>& lhs, const ArrayRef<T, Size> rhs)
{
    for (auto i=decltype(lhs.size()){0}; i<lhs.size(); ++i) {
        if (lhs[i] > rhs[i]) return false;
    }
    return true;
}

template<typename T, std::size_t Size>
bool operator>= (const ArrayRef<T, Size>& lhs, const ArrayRef<T, Size> rhs)
{
    for (auto i=decltype(lhs.size()){0}; i<lhs.size(); ++i) {
        if (lhs[i] < rhs[i]) return false;
    }
    return true;
}

template<typename T, std::size_t Size>
bool operator< (const ArrayRef<T, Size>& lhs, const ArrayRef<T, Size> rhs)
{
    return ! (lhs >= rhs);
}

template<typename T, std::size_t Size>
bool operator> (const ArrayRef<T, Size>& lhs, const ArrayRef<T, Size> rhs)
{
    return ! (lhs <= rhs);
}

template<typename T, std::size_t Index, std::size_t Size>
T& get (ArrayRef<T, Size>& array)
{
    static_assert(Index < Size, "The accessed index must be smaller than the size of the array,");
    return array[Index];
}

template<typename T, std::size_t Index, std::size_t Size>
T const& get (const ArrayRef<T, Size>& array)
{
    static_assert(Index < Size, "The accessed index must be smaller than the size of the array,");
    return array[Index];
}

template<typename T, std::size_t Size>
void swap(ArrayRef<T, Size>& a, ArrayRef<T, Size> b)
{
    a.swap(b);
}

}
