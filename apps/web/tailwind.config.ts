import type { Config } from 'tailwindcss';
import typography from '@tailwindcss/typography';

export default {
  darkMode: 'class',
  content: ['./index.html', './src/**/*.{ts,tsx}'],
  theme: {
    extend: {
      colors: {
        brand: {
          50: '#EEF2F8',
          100: '#D9E1EE',
          200: '#B3C3DD',
          300: '#7F97BF',
          400: '#4C6A9E',
          500: '#2A4B7F',
          600: '#14305C',
          700: '#10274B',
          800: '#0C1E3A',
          900: '#081528',
        },
        gold: {
          50: '#FBF7E9',
          100: '#F5EBC6',
          200: '#EBD88E',
          300: '#DFC45C',
          400: '#D4B33B',
          500: '#C9A227',
          600: '#A8861F',
          700: '#846818',
          800: '#5F4B11',
          900: '#3D300B',
        },
        steel: {
          50: '#F3F5F8',
          100: '#E4E8EE',
          200: '#C8D0DB',
          300: '#A6B2C3',
          400: '#8292A8',
          500: '#5B6B82',
          600: '#4A586C',
          700: '#3A4556',
          800: '#2B3341',
          900: '#1D232D',
        },
        surface: {
          DEFAULT: '#F5F7FA',
          card: '#FFFFFF',
          dark: '#0B1220',
          'dark-card': '#111A2E',
          'dark-border': '#223052',
        },
      },
      fontFamily: {
        sans: [
          'system-ui',
          '-apple-system',
          'Segoe UI',
          'Microsoft YaHei',
          'PingFang SC',
          'Hiragino Sans GB',
          'Noto Sans CJK SC',
          'Noto Sans SC',
          'sans-serif',
        ],
        mono: ['ui-monospace', 'SFMono-Regular', 'Menlo', 'Consolas', 'Liberation Mono', 'monospace'],
      },
      borderRadius: {
        xl: '12px',
        '2xl': '16px',
      },
      boxShadow: {
        card: '0 1px 2px rgba(16, 39, 75, 0.04), 0 4px 16px rgba(16, 39, 75, 0.06)',
        'card-hover': '0 2px 4px rgba(16, 39, 75, 0.06), 0 10px 28px rgba(16, 39, 75, 0.10)',
        'card-dark': '0 1px 2px rgba(0,0,0,0.3), 0 6px 20px rgba(0,0,0,0.35)',
      },
      keyframes: {
        'fade-in': { from: { opacity: '0', transform: 'translateY(4px)' }, to: { opacity: '1', transform: 'translateY(0)' } },
        'pulse-dot': { '0%, 100%': { opacity: '1' }, '50%': { opacity: '0.35' } },
      },
      animation: {
        'fade-in': 'fade-in 0.25s ease-out',
        'pulse-dot': 'pulse-dot 1.4s ease-in-out infinite',
      },
    },
  },
  plugins: [typography],
} satisfies Config;
