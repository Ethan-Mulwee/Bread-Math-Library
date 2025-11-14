#ifndef BML_IOSTREAM
#define BML_IOSTREAM

#include <iostream>
#include <sstream>

#include <string>
#include "bml_types.hpp"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

namespace bml {

    inline std::ostream &operator <<(std::ostream &os, const Vector2 &v) {
        os << "(" 
        ANSI_COLOR_RED + std::to_string(v.x) + ANSI_COLOR_RESET ", " 
        ANSI_COLOR_GREEN + std::to_string(v.y) + ANSI_COLOR_RESET
        ")";
        return os;
    }

    inline std::ostream &operator <<(std::ostream &os, const Vector3 &v) {
        const char* xSpace = std::signbit(v.x) ? "" : " ";
        const char* ySpace = std::signbit(v.y) ? "" : " ";
        const char* zSpace = std::signbit(v.z) ? "" : " ";
        os << "(" 
        ANSI_COLOR_RED   + (xSpace + std::to_string(v.x)) + ANSI_COLOR_RESET ", " 
        ANSI_COLOR_GREEN + (ySpace + std::to_string(v.y)) + ANSI_COLOR_RESET ", " 
        ANSI_COLOR_BLUE  + (zSpace + std::to_string(v.z)) + ANSI_COLOR_RESET
        ")";
        return os;
    }

    inline std::ostream &operator <<(std::ostream &os, const Vector4 &v) {
        os << "(" 
        ANSI_COLOR_RED + std::to_string(v.x) + ANSI_COLOR_RESET ", " 
        ANSI_COLOR_GREEN + std::to_string(v.y) + ANSI_COLOR_RESET ", " 
        ANSI_COLOR_BLUE + std::to_string(v.z) + ANSI_COLOR_RESET ", "
        ANSI_COLOR_MAGENTA + std::to_string(v.w) + ANSI_COLOR_RESET
        ")";
        return os;
    }

    inline std::ostream &operator <<(std::ostream &os, const Quaternion &q) {
        os << "(" 
        ANSI_COLOR_RED + std::to_string(q.x) + ANSI_COLOR_RESET ", " 
        ANSI_COLOR_GREEN + std::to_string(q.y) + ANSI_COLOR_RESET ", " 
        ANSI_COLOR_BLUE + std::to_string(q.z) + ANSI_COLOR_RESET ", "
        ANSI_COLOR_MAGENTA + std::to_string(q.w) + ANSI_COLOR_RESET
        "), length = " + std::to_string(q.length());
        return os;
    }

    inline std::ostream &operator <<(std::ostream &os, const Transform &t) {
        os <<
        "\n| translation: " << t.translation << "\n"
        "| rotation: " << t.rotation << "\n"
        "| scale: " << t.scale << "\n";
        return os;
    }

    inline std::ostream &operator <<(std::ostream &os, const Matrix2x2 &m) {
        unsigned longest = 0;

        std::string strings[2][2];

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                unsigned length = std::to_string(m[i][j]).length();
                if (length > longest) {
                    longest = length;
                }
            }
        }

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                std::string str = std::to_string(m[i][j]);
                while (str.length() < longest) {
                    str += "0";
                }
                strings[i][j] = str;
            }
        }

        os <<
        "\n""|" ANSI_COLOR_RED + strings[0][0] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][0] + ANSI_COLOR_RESET "| \n" +
            "|" ANSI_COLOR_RED + strings[0][1] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][1] + ANSI_COLOR_RESET "| \n";
        return os;
    }

    inline std::ostream &operator <<(std::ostream &os, const Matrix3x3 &m) {
        unsigned longest = 0;

        std::string strings[3][3];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                unsigned length = std::to_string(m[i][j]).length();
                if (length > longest) {
                    longest = length;
                }
            }
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::string str = std::to_string(m[i][j]);
                while (str.length() < longest) {
                    str += "0";
                }
                strings[i][j] = str;
            }
        }

        os <<
        "\n""|" ANSI_COLOR_RED + strings[0][0] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][0] + ANSI_COLOR_RESET ", " ANSI_COLOR_BLUE + strings[2][0] + ANSI_COLOR_RESET "| \n" +
            "|" ANSI_COLOR_RED + strings[0][1] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][1] + ANSI_COLOR_RESET ", " ANSI_COLOR_BLUE + strings[2][1] + ANSI_COLOR_RESET "| \n" +
            "|" ANSI_COLOR_RED + strings[0][2] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][2] + ANSI_COLOR_RESET ", " ANSI_COLOR_BLUE + strings[2][2] + ANSI_COLOR_RESET "| \n";
        return os;
    }

    inline std::ostream &operator <<(std::ostream &os, const Matrix4x4 &m) {
        unsigned longest = 0;

        std::string strings[4][4];

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                unsigned length = std::to_string(m[i][j]).length();
                if (length > longest) {
                    longest = length;
                }
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                std::string str = std::to_string(m[i][j]);
                while (str.length() < longest) {
                    str += "0";
                }
                strings[i][j] = str;
            }
        }

        os <<
                "\n""|" ANSI_COLOR_RED + strings[0][0] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][0] + ANSI_COLOR_RESET ", " ANSI_COLOR_BLUE + strings[2][0] + ANSI_COLOR_RESET + ", " ANSI_COLOR_MAGENTA + strings[3][0] + ANSI_COLOR_RESET "| \n" +
                    "|" ANSI_COLOR_RED + strings[0][1] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][1] + ANSI_COLOR_RESET ", " ANSI_COLOR_BLUE + strings[2][1] + ANSI_COLOR_RESET + ", " ANSI_COLOR_MAGENTA + strings[3][1] + ANSI_COLOR_RESET "| \n" +
                    "|" ANSI_COLOR_RED + strings[0][2] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][2] + ANSI_COLOR_RESET ", " ANSI_COLOR_BLUE + strings[2][2] + ANSI_COLOR_RESET + ", " ANSI_COLOR_MAGENTA + strings[3][2] + ANSI_COLOR_RESET "| \n" +
                    "|" ANSI_COLOR_RED + strings[0][3] + ANSI_COLOR_RESET ", " ANSI_COLOR_GREEN + strings[1][3] + ANSI_COLOR_RESET ", " ANSI_COLOR_BLUE + strings[2][3] + ANSI_COLOR_RESET + ", " ANSI_COLOR_MAGENTA + strings[3][3] + ANSI_COLOR_RESET "| \n";
        return os;
    }

}

#endif