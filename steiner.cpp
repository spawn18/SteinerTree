#include "SFML/Graphics.hpp"

#include <windows.h>
#include <clocale>
#include <iostream>
#include <vector>
#include <cmath>
#include <queue>
#include <array>
#include <algorithm>



constexpr double TERMINAL_RADIUS = 5.0;
constexpr double PI = 3.14159265358979323846;
constexpr double TERMINAL_LIMIT = 4;
constexpr int    SCREEN_W = 400;
constexpr int    SCREEN_H = 400;
constexpr double EPS = 1E-8;

bool is_almost_zero(double x)
{
    return fabs(x) < EPS;
}

double distance(const sf::Vector2f& p1, const sf::Vector2f& p2)
{
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

double prim(const std::vector<sf::Vector2f>& t, std::vector<sf::Vertex>& v)
{
    std::priority_queue<std::tuple<double, int, int>, std::vector<std::tuple<double, int, int>>, std::greater<std::tuple<double, int, int>>> q;


    std::vector<bool> treeContains (t.size(), false );

    for (int i = 0; i < t.size(); i++)
        q.push(std::make_tuple(distance(t[0], t[i]), 0, i));

    double cost = 0;

    while (!q.empty())
    {
        auto p = q.top(); q.pop();

        auto l = std::get<0>(p);
        auto from = std::get<1>(p);
        auto to = std::get<2>(p);

        if (!treeContains[to])
        {
            v.emplace_back(t[from], sf::Color::Black);
            v.emplace_back(t[to], sf::Color::Black);

            cost += l;
            treeContains[to] = true;

            for (int i = 0; i < t.size(); i++)
            {
                if (!treeContains[i]) q.push(std::make_tuple(distance(t[to], t[i]), to, i));
            }
        }
    }

    return cost;
}

bool is_in_triangle(const sf::Vector2f& p, const sf::Vector2f& v1, const sf::Vector2f& v2, const sf::Vector2f& v3)
{
    double d1, d2, d3;
    bool has_neg, has_pos;

    d1 = (p.x - v2.x) * (v1.y - v2.y) - (v1.x - v2.x) * (p.y - v2.y);
    d2 = (p.x - v3.x) * (v2.y - v3.y) - (v2.x - v3.x) * (p.y - v3.y);
    d3 = (p.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (p.y - v1.y);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

bool is_convex(const std::vector<sf::Vector2f>& t)
{
    return !is_in_triangle(t[0], t[1], t[2], t[3]) &&
        !is_in_triangle(t[1], t[0], t[2], t[3]) &&
        !is_in_triangle(t[2], t[0], t[1], t[3]) &&
        !is_in_triangle(t[3], t[0], t[1], t[2]);
}

bool segments_intersect(const sf::Vector2f& a, const sf::Vector2f& b, const sf::Vector2f& c, const sf::Vector2f& d)
{
    double det = (b.x - a.x) * (c.y - d.y) - (b.y - a.y) * (c.x - d.x);
    double det1 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    double det2 = (c.x - a.x) * (c.y - d.y) - (c.y - a.y) * (c.x - d.x);

    if (!is_almost_zero(det))
    {
        double t1 = det1 / det;
        double t2 = det2 / det;

        if (0 <= t1 && t1 <= 1 && 0 <= t2 && t2 <= 1)
        {
            return true;
        }
    }

    return false;
}

bool solve_quad(long double A, long double B, long double C, long double& x1, long double& x2)
{
    long double D = B * B - 4 * A * C;
    if (D >= 0)
    {
        x1 = (-B + sqrt(D)) / (2.0 * A);
        x2 = (-B - sqrt(D)) / (2.0 * A);

        return true;
    }

    return false;
}

bool arc_intersect(const sf::Vector2f& a, const sf::Vector2f& b, const sf::Vector2f& c, double rad, sf::Vector2f& f)
{
    sf::Vector2f dv = b - a;

    long double t1, t2;
    long double A = dv.x * dv.x + dv.y * dv.y;
    long double B = 2.0 * ((a.x - c.x) * dv.x + (a.y - c.y) * dv.y);
    long double C = pow(a.x - c.x, 2) + pow(a.y - c.y, 2) - rad * rad;

    bool r = solve_quad(A, B, C, t1, t2);
    if (r == true)
    {
        double s = !is_almost_zero(t1) ? t1 : t2;
        f = a + (float)s * dv;
        return true;
    }

    return false;
}

double length(const sf::Vector2f& v)
{
    return sqrt(v.x * v.x + v.y * v.y);
}

double angle_between(const sf::Vector2f& v1, const sf::Vector2f& v2)
{
    return acos((v1.x * v2.x + v1.y * v2.y) / (length(v1) * length(v2)));
}

double angle(const sf::Vector2f& v)
{
    double a = std::atan2(v.y, v.x);
    return a > 0 ? a : 2 * PI + a;
}

void rotate(sf::Vector2f& v, double angle)
{
    double x = cos(angle) * v.x - sin(angle) * v.y;
    double y = sin(angle) * v.x + cos(angle) * v.y;
    v.x = x;
    v.y = y;
}

void log(int n, double LST, double LMT)
{
    std::cout << "Кол-во: " << n;
    std::cout << " - Штейнер: " << LST;
    std::cout << " - Минимальное: " << LMT;
    std::cout << " - Отношение: " << LST / LMT;
    std::cout << " - Это на ~" << 100 * (1 - LST / LMT) << "% короче чем минимальное остовное дерево" << std::endl;
}

double tree_length(const std::vector<sf::Vertex>& tree)
{
    double L = 0;
    for (int i = 0; i < tree.size() - 1; i += 2)
        L += distance(tree[i].position, tree[i + 1].position);
    return L;
}

double steiner_n(int n, const std::vector<sf::Vector2f>& t, std::vector<sf::Vector2f>& s, std::vector<sf::Vertex>& v)
{
    if (n == 2)
    {
        v.emplace_back(t[0], sf::Color::Black);
        v.emplace_back(t[1], sf::Color::Black);

        return tree_length(v);
    }
    if (n == 3)
    {
        sf::Vector2f av = t[2] - t[1];
        sf::Vector2f bv = t[0] - t[2];
        sf::Vector2f cv = t[1] - t[0];

        double A = angle_between(cv, -bv);
        double B = angle_between(av, -cv);
        double C = angle_between(-av, bv);

        if (A < 2 * PI / 3 && B < 2 * PI / 3 && C < 2 * PI / 3)
        {
            double a = length(av);
            double b = length(bv);
            double c = length(cv);

            double x = a * (1.0 / sin(A + PI / 3.0));
            double y = b * (1.0 / sin(B + PI / 3.0));
            double z = c * (1.0 / sin(C + PI / 3.0));
            double w = a * (1.0 / sin(A + PI / 3.0)) + b * (1.0 / sin(B + PI / 3.0)) + c * (1.0 / sin(C + PI / 3.0));

            float k1 = x / w;
            float k2 = y / w;
            float k3 = z / w;

            s.emplace_back(k1 * t[0] + k2 * t[1] + k3 * t[2]);

            v.emplace_back(t[0], sf::Color::Black);
            v.emplace_back(s[0], sf::Color::Black);
            v.emplace_back(t[1], sf::Color::Black);
            v.emplace_back(s[0], sf::Color::Black);
            v.emplace_back(t[2], sf::Color::Black);
            v.emplace_back(s[0], sf::Color::Black);

            return tree_length(v);
        }
    }
    else if (n == 4)
    {
        if (is_convex(t))
        {
            sf::Vector2f v1 = t[1] - t[0];
            sf::Vector2f v2 = t[3] - t[2];

            rotate(v1, -PI / 3);
            rotate(v2, -PI / 3);

            sf::Vector2f e1 = t[0] + v1;
            sf::Vector2f e2 = t[2] + v2;

            sf::Vector2f c1((t[0].x + t[1].x + e1.x) / 3, (t[0].y + t[1].y + e1.y) / 3);
            sf::Vector2f c2((t[2].x + t[3].x + e2.x) / 3, (t[2].y + t[3].y + e2.y) / 3);

            double r1 = distance(t[0], c1);
            double r2 = distance(t[2], c2);

            if (distance(c1, c2) > r1 + r2)
            {
                if (segments_intersect(e1, e2, t[0], t[1]) && segments_intersect(e1, e2, t[2], t[3]))
                {
                    sf::Vector2f f1, f2;

                    if (arc_intersect(e1, e2, c1, r1, f1) && arc_intersect(e2, e1, c2, r2, f2))
                    {
                        s.push_back(f1);
                        s.push_back(f2);

                        v.emplace_back(t[0], sf::Color::Black);
                        v.emplace_back(s[0], sf::Color::Black);
                        v.emplace_back(t[1], sf::Color::Black);
                        v.emplace_back(s[0], sf::Color::Black);
                        v.emplace_back(s[0], sf::Color::Black);
                        v.emplace_back(s[1], sf::Color::Black);
                        v.emplace_back(t[2], sf::Color::Black);
                        v.emplace_back(s[1], sf::Color::Black);
                        v.emplace_back(t[3], sf::Color::Black);
                        v.emplace_back(s[1], sf::Color::Black);

                        return tree_length(v);
                    }
                }
            }
        }
    }

    return prim(t, v);
}

double steiner(std::vector<sf::Vector2f>& t, std::vector<sf::Vector2f>& s, std::vector<sf::Vertex>& v)
{
    if (t.size() == 2)
    {
        return steiner_n(2, t, s, v);
    }
    else if (t.size() == 3)
    {
        return steiner_n(3, t, s, v);
    }
    else if (t.size() == 4)
    {
        sf::Vector2f d = { -(t[0].x + t[1].x + t[2].x + t[3].x) / 4.f, -(t[0].y + t[1].y + t[2].y + t[3].y) / 4.f };
        std::sort(begin(t), end(t), [&d](sf::Vector2f p1, sf::Vector2f p2) { return angle(p1 + d) < angle(p2 + d); });

        std::vector<sf::Vertex> MST;
        double LMST = prim(t, MST);
        double L = LMST;

        std::vector<sf::Vector2f> steins;
        std::vector<sf::Vertex> verts;

        std::vector<std::vector<sf::Vector2f>> combinations =
        {
            {t[0],t[1],t[2],t[3]},
            {t[1],t[2],t[3],t[0]},
        };

        for (int i = 0; i < combinations.size(); i++)
        {
            verts.clear();
            steins.clear();

            double l = steiner_n(4, combinations[i], steins, verts);
            if (l < L)
            {
                v.clear(); s.clear();
                std::copy(begin(verts), end(verts), std::back_inserter(v));
                std::copy(begin(steins), end(steins), std::back_inserter(s));
                L = l;
            }
        }

        if (L < LMST) return L;

        combinations =
        {
            {t[0],t[1],t[2]},
            {t[1],t[2],t[3]},
            {t[2],t[3],t[0]},
            {t[3],t[0],t[1]}
        };

        L = LMST;

        for (int i = 0; i < combinations.size(); i++)
        {
            verts.clear();
            steins.clear();

            steiner_n(3, combinations[i], steins, verts);

            sf::Vector2f p1, p2;

            for (int j = 0; j < t.size(); j++)
            {
                auto res = std::find(begin(combinations[i]), end(combinations[i]), t[j]);
                if (res == end(combinations[i]))
                {
                    p1 = t[j];
                    break;
                }
            }

            std::vector<std::pair<double, sf::Vector2f>> dists;

            for (int j = 0; j < combinations[i].size(); j++)
                dists.emplace_back(distance(p1, combinations[i][j]), combinations[i][j]);

            std::sort(begin(dists), end(dists), [](const auto& p1, const auto& p2) { return p1.first < p2.first; });
            p2 = dists[0].second;

            verts.emplace_back(p1, sf::Color::Black);
            verts.emplace_back(p2, sf::Color::Black);

            double l = tree_length(verts);
            if (l < L)
            {
                s.clear(); v.clear();
                std::copy(begin(steins), end(steins), std::back_inserter(s));
                std::copy(begin(verts), end(verts), std::back_inserter(v));
                L = l;
            }
        }

        if (L == LMST)
            std::copy(begin(MST), end(MST), std::back_inserter(v));

        return L;
    }
}

int main()
{
    setlocale(LC_ALL, "Russian");
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

    std::string title = u8"Минимальное евклидовое дерево Штейнера";
    sf::RenderWindow window(sf::VideoMode(400, 400), sf::String::fromUtf8(title.begin(), title.end()), sf::Style::Titlebar | sf::Style::Close);

    std::vector<sf::Vector2f> terminals;
    std::vector<sf::Vector2f> steiners;
    std::vector<sf::Vertex> SMT, MST;

    bool update = false;

    sf::CircleShape t;
    t.setOrigin(TERMINAL_RADIUS, TERMINAL_RADIUS);
    t.setFillColor(sf::Color::Blue);
    t.setRadius(TERMINAL_RADIUS);

    sf::CircleShape s;
    s.setOrigin(TERMINAL_RADIUS / 2, TERMINAL_RADIUS / 2);
    s.setFillColor(sf::Color::Black);
    s.setRadius(TERMINAL_RADIUS / 2);

    sf::Clock clock;
    int flash = 1;

    while (window.isOpen())
    {
        sf::Event e;
        while (window.pollEvent(e))
        {
            switch (e.type)
            {
            case sf::Event::Closed:
            {
                window.close();
                break;
            }

            case sf::Event::MouseButtonPressed:
            {
                if (e.mouseButton.button == sf::Mouse::Left)
                {
                    if (terminals.size() < TERMINAL_LIMIT)
                    {
                        sf::Vector2f point(e.mouseButton.x, e.mouseButton.y);
                        if (find(begin(terminals), end(terminals), point) == end(terminals))
                        {
                            terminals.push_back(point);
                            if (terminals.size() >= 2)
                                update = true;
                        }
                    }
                }
                break;
            }
            }
        }

        window.clear(sf::Color::White);

        if (update)
        {
            SMT.clear();
            MST.clear();
            steiners.clear();

            double lenSMT = steiner(terminals, steiners, SMT);
            double lenMST = prim(terminals, MST);
            std::for_each(begin(MST), end(MST), [](sf::Vertex& v) { v.color = sf::Color::Green; });

            log(terminals.size(), lenSMT, lenMST);
            update = false;
        }

        window.draw(SMT.data(), SMT.size(), sf::Lines);

        if (clock.getElapsedTime().asMilliseconds() >= 750)
        {
            flash = 1 - flash;
            clock.restart();
        }

        if (flash) window.draw(MST.data(), MST.size(), sf::Lines);

        for (const auto& i : steiners)
        {
            s.setPosition(i.x, i.y);
            window.draw(s);
        }

        for (const auto& i : terminals)
        {
            t.setPosition(i.x, i.y);
            window.draw(t);
        }

        window.display();
    }

    return 0;
}
