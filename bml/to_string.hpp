#ifndef BML_STRING
#define BML_STRING

#include <string>
#include "smath_types.hpp"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

namespace bml {

    std::string to_string(const Vector2 &v) {
        return "(" + std::to_string(v.x) + "," + std::to_string(v.y) + ")";
    }

    std::string to_string(const Vector3 &v) {
        return "(" + std::to_string(v.x) + "," + std::to_string(v.y) + "," + std::to_string(v.z) + ")";
    }

    std::string to_string(const Vector4 &v) {
        return "(" + std::to_string(v.x) + "," + std::to_string(v.y) + "," + std::to_string(v.z) + "," + std::to_string(v.w) + ")";
    }

    std::string to_string(const Matrix2x2 &m) {
        return
            "|" + std::to_string(m[0][0]) + ", " + std::to_string(m[0][1]) + "| \n" +
            "|" + std::to_string(m[1][0]) + ", " + std::to_string(m[1][1]) + "| \n";
    }

    std::string to_string(const Matrix3x3 &m) {
        return 
            "|" + std::to_string(m[0][0]) + ", " + std::to_string(m[0][1]) + ", " + std::to_string(m[0][2]) + "| \n" +
            "|" + std::to_string(m[1][0]) + ", " + std::to_string(m[1][1]) + ", " + std::to_string(m[1][2]) + "| \n" +
            "|" + std::to_string(m[2][0]) + ", " + std::to_string(m[2][1]) + ", " + std::to_string(m[2][2]) + "| \n";
    }

    std::string to_string(const Matrix4x4 &m) {
        return
            "|" + std::to_string(m[0][0]) + ", " + std::to_string(m[0][1]) + ", " + std::to_string(m[0][2]) + ", " + std::to_string(m[0][3]) + "| \n" +
            "|" + std::to_string(m[1][0]) + ", " + std::to_string(m[1][1]) + ", " + std::to_string(m[1][2]) + ", " + std::to_string(m[1][3]) + "| \n" +
            "|" + std::to_string(m[2][0]) + ", " + std::to_string(m[2][1]) + ", " + std::to_string(m[2][2]) + ", " + std::to_string(m[2][3]) + "| \n" +
            "|" + std::to_string(m[3][0]) + ", " + std::to_string(m[3][1]) + ", " + std::to_string(m[3][2]) + ", " + std::to_string(m[3][3]) + "|";
    }

    std::string to_string_pretty(const Matrix3x3 &m) {
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

        return
            "|" + strings[0][0] + ", " + strings[0][1] + ", " + strings[0][2] + "| \n" +
            "|" + strings[1][0] + ", " + strings[1][1] + ", " + strings[1][2] + "| \n" +
            "|" + strings[2][0] + ", " + strings[2][1] + ", " + strings[2][2] + "| \n";
    }

    std::string to_string_pretty(const Matrix4x4 &m) {
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

        return
            "|" + strings[0][0] + ", " + strings[0][1] + ", " + strings[0][2] + ", " + strings[0][3] + "| \n" +
            "|" + strings[1][0] + ", " + strings[1][1] + ", " + strings[1][2] + ", " + strings[1][3] + "| \n" +
            "|" + strings[2][0] + ", " + strings[2][1] + ", " + strings[2][2] + ", " + strings[2][3] + "| \n" +
            "|" + strings[3][0] + ", " + strings[3][1] + ", " + strings[3][2] + ", " + strings[3][3] + "|";
    }

    std::string to_string(const Quaternion &q) {
        return "( x:" + std::to_string(q.x) + ", y:" + std::to_string(q.y) + ", z:" + std::to_string(q.z) + ", w:" + std::to_string(q.w) + ")";
    }

    std::string to_string_verbose (const Quaternion &q) {
        return std::string(
            "Quaternion: (" ANSI_COLOR_RED "x:" + std::to_string(q.x) + ANSI_COLOR_RESET "," ANSI_COLOR_GREEN " y:" + std::to_string(q.y) + ANSI_COLOR_RESET "," ANSI_COLOR_BLUE " z:" + std::to_string(q.z) + ANSI_COLOR_RESET "," ANSI_COLOR_YELLOW " w:" + std::to_string(q.w) + ANSI_COLOR_RESET ") " + 
            "Length: (" + std::to_string(q.length()) + "), Axis: " + to_string(q.axis()) + ", Angle: (" + std::to_string(q.angle()) + ")"
        );
    }
}

#endif