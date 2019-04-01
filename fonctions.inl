
template <typename T>
void shuffle_n(std::vector<T>& vec, int n) {
    int i = 0;
    while (n--) {
        std::swap(vec[i], vec[randomInt(i, vec.size()-1)]);
        i++;
    }
}

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> vec) {
    out << '{';
    for (const T& elem : vec) {
        out << elem << ',';
    }
    out << '}';
    return out;
}